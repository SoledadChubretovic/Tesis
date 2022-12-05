#### NSGA II ####

import numpy as np
import time
from Thermal_1 import ThermalAnalysis
from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.rnd import IntegerRandomSampling
from pymoo.operators.repair.rounding import RoundingRepair
from pymoo.algorithms.moo.nsga2 import RankAndCrowdingSurvival
from pymoo.optimize import minimize
from tabulate import tabulate
import matplotlib.pyplot as plt
from pymoo.termination.default import DefaultMultiObjectiveTermination
#from pymoo.termination import get_termination

from functions import (
    generate_cavitiy_ansys_parameters,
    tabique_w,
    tabique_l,
    upper_limit_l
)
from k_constants import (
    N_CAVITIES,
    N_VARIABLES,
    L,
    H,
    T_MIN,
    W2_MIN,
    W2_MAX,
    LAMBDA_CLAY100_MIN,
    LAMBDA_CLAY100_MAX,
    w_MAX,
    N_CONSTRAINTS,
    N_OBJECTIVES,
    G2,
    #N_GENERATION,
    N_OFFSPRINGS,
    POPULATION_SIZE,
    CAVITIES_PER_ROW,
    U_MURO_INVALID
)

start = time.time()

f = open(r"C:\Users\nchubretovic\OneDrive - Entel\Escritorio\Sole\Tesis\1_Model\Error.txt", "w")
f.close()

#### IMPLEMENT THE PROBLEM ####

# number of variables, constrains and objective functions
n_var = N_VARIABLES #number of variables (li,w,W,lambda_clay)
n_ieq_constr = N_CONSTRAINTS #number of constraints
n_obj = N_OBJECTIVES #number of objectives


# lower limit fo variables
xl = np.ones(N_VARIABLES - 2) * T_MIN  # mm
xl = np.append(xl, W2_MIN)  # mm
xl = np.append(xl, LAMBDA_CLAY100_MIN)  # mm

# upper limits of variables
xu =[]
for i in range(0, len(CAVITIES_PER_ROW)):
    l_max = upper_limit_l(i)
    xu_l = np.ones(CAVITIES_PER_ROW[i]) * l_max
    xu = np.append(xu, xu_l)
xu = np.append(xu, w_MAX)  # mm
xu = np.append(xu, W2_MAX)  # mm
xu = np.append(xu, LAMBDA_CLAY100_MAX)  # mm

# ElementwiseProblem evaluates one solution at a time
class ConstrainedProblem(ElementwiseProblem):

    def __init__(self, **kwargs):
        super().__init__(n_var = n_var,
                         n_obj = n_obj,
                         n_ieq_constr = n_ieq_constr,
                         xl = xl,
                         xu = xu,
                         **kwargs)

    def _evaluate(self, x, out, *args, **kwargs):
        # geometrical parameters
        # vector x is defined as [l1, l2, ..., ln-1, ln, w, W, λ]

        _w_pos = len(x) - 3
        _W_pos = len(x) - 2
        _λ_pos = len(x) - 1

        # width of the brick
        W = x[_W_pos] * 2  # mm

        # width of every cavity
        w = x[_w_pos]  # mm

        # sum of cavitie areas within the brick
        # void area (area de huecos in spanish)
        Ah = x[_w_pos] * sum(x[:N_CAVITIES]) #mm2

        # total area of the brick
        # gross area (area bruta in spanish)
        Ab = W * L #mm2

        # net area of the brick
        An = Ab - Ah #mm2

        # thermal conductivity of the clay
        λ_clay = x[_λ_pos] / 100  # W/mK

        # structural parameters
        # density of the clay in kg/mm3
        ρ_clay = 2.0841e-6*λ_clay + 0.55042e-6 # kg/mm3
        # compression resistance of the clay in MPa
        compression_resistance_clay = 76.766*λ_clay-26.362 # MPa
        # compression resistance of the brick in MPa
        compression_resistance_brick = compression_resistance_clay*An/Ab # MPa
        # compression resistance limit of the brick in MPa
        # based in NCh1928
        compression_resistance_limit = G2 #MPa

        # Calculo U_muro
        # if partition walls are smaller than the limit U = 50 (big number)
        # if net area is less than 50% of gross area U = 50 (big number)
        # this is done to reduce running time to
        # avoid using ANSYS when the geometry  of the brick is wrong
        # or it doesnt fulfill the constraints
        U_muro = 0

        thickness_tabique_w = tabique_w(W,w)
        if thickness_tabique_w < T_MIN:
            U_muro = U_MURO_INVALID
        
        for i in range(0, len(CAVITIES_PER_ROW)):
            thickness_tabique_l = tabique_l(i, x)
            if thickness_tabique_l < T_MIN:
                U_muro = U_MURO_INVALID
                break

        if An/Ab < 0.5:
            U_muro = U_MURO_INVALID
        
        if U_muro == 0:
            # needed to create brick in ANSYS
            # (center_x, center_y, dimension_x, dimension_y) of perimeter of the brick
            b = np.array([L/2, W/2, L, W])/1000  # m

            # (center_x, center_y, dimension_x, dimension_y) of each cavity
            matrix = generate_cavitiy_ansys_parameters(x, W, w) / 1000  # m
            U_muro = ThermalAnalysis(b, matrix, λ_clay)

        # Objective 1: minimize weigth
        solid_volume = H * An # mm3
        f1 = ρ_clay * solid_volume # kg

        # Objective 2: minimize wall heat transfer
        f2 = U_muro # W/mK

        # Constraint 1: structural constraint
        # compression resistance of the brick must be higher than the limit
        g1 = - compression_resistance_brick + compression_resistance_limit

        # Constraint 2: void area must be less than 50% of gross area
        g2 = Ah/Ab - 0.5

        # Constraint 3: must be a valid geometry
        # up in the code when partition wall thickness is less than 10 mm
        # U_muro is set to 50
        # we can use that value to set the geometrical constraint
        g3 = U_muro - (U_MURO_INVALID - 1)

        out["F"] = [f1,f2]
        out["G"] = [g1,g2,g3]

problem = ConstrainedProblem()

#### INITIALIZE ALGORITHM NSGA2 ####

# algorithm parameters based in ZDT benchmark problems
algorithm = NSGA2(
    pop_size = POPULATION_SIZE, # population size
    n_offsprings = N_OFFSPRINGS, # new individuals from crosssover and mutation (children),
    sampling = IntegerRandomSampling(), # creates initial population
    crossover = SBX(prob = 0.9, eta = 15, vtype = float, repair = RoundingRepair()),
    mutation = PM(prob = 1 / n_var, eta = 20, vtype = float, repair = RoundingRepair()),
    survival = RankAndCrowdingSurvival(),
    eliminate_duplicates = True, # element created are different from the ones that already exist
)

#### DEFINE A TERMINATION CRITERION ####

# DefaultMultiObjectiveTermination has the following parameters:
#     xtol = 0.0005,
#     cvtol = 1e-8,
#     ftol = 0.005,
#     period = 50,
#     n_max_gen = 1000,
#     n_max_evals = 100000
termination = DefaultMultiObjectiveTermination()

class MyMultiObjectiveDefaultTermination(DefaultMultiObjectiveTermination):

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.reasons = None

    def _decide(self, metrics):
        decisions = metrics[-1]
        cont = decisions["x_tol"] and (decisions["cv_tol"] or decisions["f_tol"])

        if not cont:
            self.reasons = [name for name in ["x_tol", "cv_tol", "f_tol"] if not decisions[name]]

        return cont

ret = minimize(problem,
               algorithm,
               termination=MyMultiObjectiveDefaultTermination(),
               seed=1,
               save_history=False,
               verbose=True)

print(ret.algorithm.termination.reasons)

#### OPTIMIZE ####

res = minimize(problem,
               algorithm,
               termination,
               seed = 1, # ramdom parameter to ensure reproducible results
               save_history = True,
               verbose = True # some printouts during the algorithm´s execution is provided
)


#### WRITE RESULTS IN TXT FILE ####
# res.F is a list size NxO with the objective values for each optimal solution
# res.X is a list size NxV with the variable values for each optimal solution
# N is the number of optimal solutions of the pareto front
# O is the number of objectives
# V is the number of variables
data1 = res.F
data2 = res.X

table_data = []
for i in range(0,len(data1)):
    linedata = []
    for j in data1[i]:
        linedata.append(j)
    for k in data2[i]:
        linedata.append(k)
    linedata[-1] = linedata[-1]/100
    linedata[-2] = linedata[-2]*2
    table_data.append(linedata)

# column tags
# column_names = ["Weight [kg]","U_value [W/m2K]"]
# for i in range(0,N_VARIABLES - 3):
#     column_names = column_names + (["l" + str(i + 1) + " [mm]"])
# column_names = column_names + (["w [mm]"])
# column_names = column_names + (["W [mm]"])
# column_names = column_names + (["λ_clay [W/mK]"])

result_table = tabulate(table_data, tablefmt = "plain")

# write table to results file
results = open(r"C:\Users\nchubretovic\OneDrive - Entel\Escritorio\Sole\Tesis\1_Model\OptPython_Log.txt", "w")
results.write(result_table)
results.close()

end = time.time()
print("Time: ", (end-start), "s")

#### VISUALIZE ####

# plot optimal pareto front
plt.figure(figsize=(7, 5))
plt.scatter(res.F[:,0], res.F[:,1], s=30, facecolors='none', edgecolors='r')
plt.xlabel("Weight kg")
plt.ylabel("Thermal Transmitance W/m2K")
plt.title("Weight vs. Thermal Transmitance")
plt.show()
plt.savefig(r"C:\Users\nchubretovic\OneDrive - Entel\Escritorio\Sole\Tesis\1_Model\Pareto_images\Pareto_" + time.strftime("%Y-%m-%d %H.%M.%S") + ".png")
