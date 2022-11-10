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
from pymoo.termination import get_termination
from pymoo.optimize import minimize
from tabulate import tabulate
import matplotlib.pyplot as plt

from functions import (
    generate_cavitiy_ansys_parameters
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
    l_MAX,
    w_MAX,
    N_CONSTRAINTS,
    N_OBJECTIVES,
    G2,
    N_GENERATION,
    N_OFFSPRINGS,
    POPULATION_SIZE
)

start = time.time()

f = open(r"C:\Users\nchubretovic\OneDrive - Entel\Escritorio\Sole\Tesis\OptPython\1_Model\Error.txt", "w")
f.close()

#### IMPLEMENT THE PROBLEM ####
    
# Number of variables, constrains and objective functions
n_var = N_VARIABLES #number of variables (li,w,W,lambda_clay)
n_ieq_constr = N_CONSTRAINTS #number of constraints
n_obj = N_OBJECTIVES #number of objectives


# lower limit fo variables
xl = np.ones(N_VARIABLES - 2) * T_MIN  # mm
xl = np.append(xl, W2_MIN)  # mm
xl = np.append(xl, LAMBDA_CLAY100_MIN)  # mm

# upper limits of variables
xu = np.ones(N_CAVITIES) * l_MAX  # mm
xu = np.append(xu, w_MAX)  # mm
xu = np.append(xu, W2_MAX)  # mm
xu = np.append(xu, LAMBDA_CLAY100_MAX)  # mm

class ConstrainedProblem(ElementwiseProblem): #ElementwiseProblem evaluates one solution at a time
    
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

        # area of the cavities within the brick ()
        Ah = x[_w_pos] * sum(x[:N_CAVITIES])

        # total area of the brick
        Ab = W * L

        # net area of the brick (total - cavities)
        An = Ab - Ah

        # needed to create brick in ANSYS
        # (center_x, center_y, dimension_x, dimension_y) of perimeter of the brick
        b = np.array([L/2, W/2, L, W])/1000  # m

        # (center_x, center_y, dimension_x, dimension_y) of each cavity
        matrix = generate_cavitiy_ansys_parameters(x, W, w) / 1000  # m

        # thermal conductivity of the clay
        λ_clay = x[_λ_pos] / 100  # W/mK

        # Structural parameters
        ρ_clay = 2.0841e-6*λ_clay + 0.55042e-6 #kg/mm3
        compression_resistance_clay = 76.766*λ_clay-26.362 #MPa
        compression_resistance_brick = compression_resistance_clay*An/Ab #MPa
        compression_resistance_limit = G2 #Mpa

        # Calculo U_muro
        if tla < 10 or tlb < 10 or tlc < 10 or tld < 10 or tw < 10:
            U_muro = 50
            if Ah/Ab < 0.5:
                U_muro = 50
            else: 
                U_muro = ThermalAnalysis(x)
                f = open(r"C:\Users\nchubretovic\OneDrive - Entel\Escritorio\Sole\Tesis\OptPython\1_Model\Error.txt", "a")
                f.write(U_muro + "/n")
                f.close

        else:
            U_muro = ThermalAnalysis(x)

        # Objective 1: minimize weigth
        f1 = ρ_clay*H*An #kg

        # Objective 2: minimize wall heat transfer
        f2 = U_muro #W/mK

        # Structural Constraint
        g1 = - compression_resistance_brick + compression_resistance_limit

        # Geometrical Constraint
        g2 = 10 - tla
        g3 = 10 - tlb
        g4 = 10 - tlc
        g5 = 10 - tld
        g6 = 10 - tw

        g7 = Ah/Ab - 0.5 # %huecos menor o igual a 50%

        out["F"] = [f1,f2]
        out["G"] = [g1,g2,g3,g4,g5,g6,g7]

problem = ConstrainedProblem()

#### INITIALIZE ALGORITHM NSGA2 ####

# PARAMETROS BSADOS EN ZDT
algorithm = NSGA2(
    pop_size = POPULATION_SIZE, #population size
    n_offsprings = N_OFFSPRINGS, #new individuals from crosssover and mutation (children),
    sampling = IntegerRandomSampling(), #creates initial population
    crossover = SBX(prob=0.9, eta=15, vtype=float, repair=RoundingRepair()),
    mutation = PM(prob=1/n_var, eta=20, vtype=float, repair=RoundingRepair()),
    survival = RankAndCrowdingSurvival(), #Yo Agregue esto
    eliminate_duplicates = True, #element created are different from the ones that already exist
)

#### DEFINE A TERMINATION CRITERION ####

termination = get_termination("n_gen",N_GENERATION)

#### OPTIMIZE ####

res = minimize(problem,
               algorithm,
               termination,
               seed = 1, #ramdom parameter to ensure reproducible results
               save_history = True,
               verbose = True #some printouts during the algorithm´s execution is provided
)


#### WRITE RESULTS IN TXT FILE ####

data1 = res.F
data2 = res.X

table_data = []
for i in range(0,len(data1)):
    linedata = []
    for j in data1[i]:
        linedata.append(j)
    for k in data2[i]:
        linedata.append(k)
    table_data.append(linedata)

#col_names = ["Weight [kg]","U_value [W/m2K]","l1 [mm]","l2 [mm]","l3 [mm]","l4 [mm]","w [mm]", "W [mm]"]
tabla = tabulate(table_data)

results = open(r"C:\Users\nchubretovic\OneDrive - Entel\Escritorio\Sole\Tesis\OptPython\1_Model\OptPython_Log.txt", "w")
results.write(tabla)
results.close()

end = time.time()
print("Time: ", (end-start), "s")

#### VISUALIZE ####

plt.figure(figsize=(7, 5))
plt.scatter(res.F[:,0], res.F[:,1], s=30, facecolors='none', edgecolors='r')
plt.xlabel("Weight kg")
plt.ylabel("Thermal Transmitance W/m2K")
plt.title("Weight vs. Thermal Transmitance")
plt.show()



