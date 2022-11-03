#### NSGA II ####

import numpy as np
from tabulate import tabulate
from Termico2 import ThermalAnalysis

#### IMPLEMENT THE PROBLEM ####

from pymoo.core.problem import ElementwiseProblem
results = open(r"C:\Users\nchubretovic\OneDrive - Entel\Escritorio\Sole\OptPython\OptPython_Log.txt", "w")
f = open(r"C:\Users\nchubretovic\OneDrive - Entel\Escritorio\Sole\OptPython\Error.txt", "w")
f.close()
# Number of variables, constrains and objective functions
n_var = 8
n_ieq_constr = 4
n_obj=2

# Upper and lower limit
xl = np.ones(n_var)*10 #mm
xu = np.array([160,160,160,160,60,60,60,60]) #mm

class ConstrainedProblem(ElementwiseProblem): #ElementwiseProblem evaluates one solution at a time
    
    def __init__(self, **kwargs):
        super().__init__(n_var = n_var,
                         n_obj = n_obj,
                         n_ieq_constr = n_ieq_constr,
                         xl = xl,
                         xu = xu,
                         **kwargs)

    def _evaluate(self, x, out, *args, **kwargs):

        #### GEOMETRY AND MATERIAL VALUES ####

        l1 = x[0]
        l2 = x[1]
        l3 = x[2]
        l4 = x[3]
        w1 = x[4]
        w2 = x[5]
        w3 = x[6]
        w4 = x[7]

        W = 100 #mm
        L = 200 #mm
        H = 20 #mm
        t_min = 10 #mm ancho minimo de tabiques

        # Areas
        Ah = w1*l1+w2*l2+w3*l3+w4*l4
        Ab = W*L
        An = abs(Ab-Ah) #mm2 porque en casos de que Ahuecos mayor que Abruta An era negativo

        # Asumo que todos los tabiques son igual de ancho en una misma linea
        tl12 = (L-l1-l2)/3 #mm
        tl34 = (L-l3-l4)/3 #mm
        tw13 = (W-w1-w3)/3 #mm
        tw24 = (W-w2-w4)/3 #mm

        # Structural parameters
        rho_clay = 1.6516e-6 #kg/mm3
        p_clay = 30 #MPa
        p_brick = p_clay*An/Ab #MPa
        p_brickG2 = 10 #Mpa

        # Thermal parameters
        k = 0.5 #W/mK
        U_Santiago = 5 #W/m2K

        # [centro x, centro y, largo x, largo y]
        b = np.array([L/2,W/2,L,W])/1000 #m
        c1 = np.array([tl12+l1/2,tw13+w1/2,l1,w1])/1000 #m
        c2 = np.array([2*tl12+l1+l2/2,tw24+w2/2,l2,w2])/1000 #m
        c3 = np.array([tl34+l3/2,2*tw13+w1+w3/2,l3,w3])/1000 #m
        c4 = np.array([2*tl34+l3+l4/2,2*tw24+w2+w4/2,l4,w4])/1000 #m
        tmin = t_min/1000 #m
        dimz = H/1000 #m

        # Calculo U_muro
        U_muro = ThermalAnalysis(tmin,k,b,c1, c2,c3,c4,dimz)

        # Objective 1: minimize weigth
        f1 = rho_clay*H*An #kg

        # Objective 3: minimize wall heat transfer
        f2 = U_muro #W/mK
        #print(U_muro)

        # Structural Constraint
        g1 = - p_brick + p_brickG2

        # Geometrical Constraints (que esten dentro del ladrillo y que no se traslapen)
        if (c1[2]+c2[2]+3*tmin) - b[2] > 0 or (c3[2]+c4[2]+3*tmin) - b[2] > 0 or (c1[3]+c3[3]+3*tmin) - b[3] > 0 or (c2[3]+c4[3]+3*tmin) - b[3] > 0:
            GP2 = 1
        else:
            GP2 = -1
        g2 = GP2

        if tl12+t_min+l1 > 2*tl34+l3 and tw13+w1+t_min > 2*tw24+w2:
            GP3 = 1
        elif tl34+t_min+l3 > 2*tl12+l1 and tw13+w3+t_min > 2*tw24+w4:
            GP3 = 1
        else:
            GP3 = -1
        g3 = GP3

        # Thermal Constrain
        g4 = U_muro - U_Santiago

        out["F"] = [f1,f2]
        out["G"] = [g1,g2,g3,g4]

problem = ConstrainedProblem()

#### INITIALIZE ALGORITHM NSGA2 ####

from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.rnd import IntegerRandomSampling
from pymoo.operators.repair.rounding import RoundingRepair
from pymoo.operators.selection.tournament import TournamentSelection
from pymoo.algorithms.moo.nsga2 import RankAndCrowdingSurvival
from pymoo.util.display.multi import MultiObjectiveOutput


# PARAMETROS BSADOS EN ZDT
algorithm = NSGA2(
    pop_size = 50, #population size
    n_offsprings = 50, #new individuals from crosssover and mutation (children),
    sampling = IntegerRandomSampling(), #creates initial population
    crossover = SBX(prob=1, eta=30, vtype=float, repair=RoundingRepair()),
    mutation = PM(prob=1/20, eta=20, vtype=float, repair=RoundingRepair()), #YO AGREGUE PROB prob = 1/20, 
    survival = RankAndCrowdingSurvival(), #Yo Agregue esto
    eliminate_duplicates = True, #element created are different from the ones that already exist
)


#### DEFINE A TERMINATIN CRITERION ####

# from pymoo.termination.default import DefaultMultiObjectiveTermination

# termination = DefaultMultiObjectiveTermination(
#     xtol = 1e-8, #design space tolerance
#     cvtol = 1e-6,
#     ftol = 0.0025, #objective space tolerance
#     period = 30,
#     n_max_gen = 50,
#     n_max_evals = 1000
# )

from pymoo.termination import get_termination

termination = get_termination("n_gen",70)


#### OPTIMIZE ####
from pymoo.optimize import minimize

res = minimize(problem,
               algorithm,
               termination,
               seed=1, #ramdom parameter to ensure reproducible results
               save_history=True,
               verbose=True #some printouts during the algorithm´s execution is provided
)

# print(res.X)
# print(res.F)
# print(res.G)

#### WRITE RESULTS IN TXT FILE ####

from tabulate import tabulate

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

col_names = ["Weight [kg]","U_value [W/m2K]","l1 [mm]","l2 [mm]","l3 [mm]","l4 [mm]","w1 [mm]","w2 [mm]","w3 [mm]","w4 [mm]"]
tabla = tabulate(table_data,headers=col_names)

results.write(tabla)
results.close()

#### VISUALIZE ####

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# # 3D
# fig = plt.figure(figsize=(7, 5))
# ax = fig.add_subplot(111, projection = "3d")
# ax.scatter3D(F[:,0], F[:,1],F[:,2])
# ax.set_xlabel("Peso kg")
# ax.set_ylabel("Resistencia a compresion MPa")
# ax.set_zlabel("Transmitancia W/mK")
# plt.title("3D")
# plt.show()

# iterarion_number = list(range(1, len_F+1))

#Peso y termico
plt.figure(figsize=(7, 5))
plt.scatter(res.F[:,0], res.F[:,1], s=30, facecolors='none', edgecolors='r')
plt.xlabel("Peso kg")
plt.ylabel("Transmitancia W/mK")
plt.title("peso y transmitancia")
plt.show()

# #Peso y compresion
# fig = plt.figure(figsize=(7, 5))
# plt.scatter(F[:, 0], F[:, 1], s=30, facecolors='none', edgecolors='r')
# plt.xlabel("Peso kg")
# plt.ylabel("Resistencia a compresion MPa")
# plt.title("peso y resistencia a compresion")
# plt.show()



