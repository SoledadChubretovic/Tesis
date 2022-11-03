#### NSGA II ####

import numpy as np
import time
from Therm_9_4 import ThermalAnalysis
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

start = time.time()

results = open(r"C:\Users\nchubretovic\OneDrive - Entel\Escritorio\Sole\OptPython\OptPython_Log.txt", "w")
f = open(r"C:\Users\nchubretovic\OneDrive - Entel\Escritorio\Sole\OptPython\Error.txt", "w")
f.close()

#### IMPLEMENT THE PROBLEM ####

# Number of variables, constrains and objective functions
n_cav_L = 9
n_cav_W = 4
n_cav = n_cav_L*n_cav_W
n_var = n_cav + 3 #number of variables (li,w,W,lambda_clay)
n_ieq_constr = 7 #number of constraints
n_obj = 2 #number of objectives

# Upper and lower limit
L = 320 #mm
W2_min = 140/2 #mm 70*2 = 140
W2_max = 300/2 #mm 150*2 = 300
t_min = 10 #mm
l_max = (L - 2*n_cav_L*t_min) #mm
w_max = (W2_max - 2*n_cav_W*t_min) #mm
lambda_clay100_min = 50 #W/mK
lambda_clay100_max = 62 #W/mK
xl = np.ones(n_var - 2)*t_min #mm
xl = np.append(xl,W2_min) #mm
xl = np.append(xl,lambda_clay100_min) #mm
xu = np.ones(n_cav)*l_max #mm
xu = np.append(xu,w_max) #mm
xu = np.append(xu,W2_max) #mm
xu = np.append(xu,lambda_clay100_max) #mm

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
        L = 320 #mm (Hacer calzar con L)
        H = 113 #mm
        t_min = 10 #mm ancho minimo de tabiques (Hacer calzar con xl)
        W = x[-2]*2 #mm

        # Areas
        Ah = x[-3]*(sum(x)-x[-1]-x[-2]-x[-3])
        Ab = W*L
        
        # Asumo que todos los tabiques son igual de ancho en una misma linea
        suma = 0
        for i in range(0, n_cav_L):    
            suma = suma + x[i]
        tla = (L-suma)/(n_cav_L+1) #mm

        suma = 0
        for i in range(n_cav_L, n_cav_L*2):    
            suma = suma + x[i]
        tlb = (L-suma)/(n_cav_L+1) #mm

        suma = 0
        for i in range(n_cav_L*2, n_cav_L*3):    
            suma = suma + x[i]
        tlc = (L-suma)/(n_cav_L+1) #mm

        suma = 0
        for i in range(n_cav_L*3, n_cav_L*4):    
            suma = suma + x[i]
        tld = (L-suma)/(n_cav_L+1) #mm

        tw = (W-n_cav_W*x[-3])/(n_cav_W+1) #mm
        
        # [centro x, centro y, largo x, largo y]
        b = np.array([L/2,W/2,L,W])/1000 #m
        dimz = H/1000 #m

        c1a = np.array([tla+x[0]/2 , W-(tw+x[-3]/2) , x[0] , x[-3]])/1000 #m
        c2a = np.array([2*tla+x[0]+x[1]/2 , W-(tw+x[-3]/2) , x[1] , x[-3]])/1000 #m
        c3a = np.array([3*tla+x[0]+x[1]+x[2]/2 , W-(tw+x[-3]/2) , x[2] , x[-3]])/1000 #m
        c4a = np.array([4*tla+x[0]+x[1]+x[2]+x[3]/2 , W-(tw+x[-3]/2) , x[3] , x[-3]])/1000 #m
        c5a = np.array([5*tla+x[0]+x[1]+x[2]+x[3]+x[4]/2 , W-(tw+x[-3]/2) , x[4] , x[-3]])/1000 #m
        c6a = np.array([6*tla+x[0]+x[1]+x[2]+x[3]+x[4]+x[5]/2 , W-(tw+x[-3]/2) , x[5] , x[-3]])/1000 #m
        c7a = np.array([7*tla+x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]/2 , W-(tw+x[-3]/2) , x[6] , x[-3]])/1000 #m
        c8a = np.array([8*tla+x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]/2 , W-(tw+x[-3]/2) , x[7] , x[-3]])/1000 #m
        c9a = np.array([9*tla+x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]/2 , W-(tw+x[-3]/2) , x[8] , x[-3]])/1000 #m

        c1b = np.array([tlb+x[0]/2 , W-(2*tw+x[-3]*3/2) , x[9] , x[-3]])/1000 #m
        c2b = np.array([2*tlb+x[0]+x[1]/2 , W-(2*tw+x[-3]*3/2) , x[10] , x[-3]])/1000 #m
        c3b = np.array([3*tlb+x[0]+x[1]+x[2]/2 , W-(2*tw+x[-3]*3/2) , x[11] , x[-3]])/1000 #m
        c4b = np.array([4*tlb+x[0]+x[1]+x[2]+x[3]/2 , W-(2*tw+x[-3]*3/2) , x[12] , x[-3]])/1000 #m
        c5b = np.array([5*tlb+x[0]+x[1]+x[2]+x[3]+x[4]/2 , W-(2*tw+x[-3]*3/2) , x[13] , x[-3]])/1000 #m
        c6b = np.array([6*tlb+x[0]+x[1]+x[2]+x[3]+x[4]+x[5]/2 , W-(2*tw+x[-3]*3/2) , x[14] , x[-3]])/1000 #m
        c7b = np.array([7*tlb+x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]/2 , W-(2*tw+x[-3]*3/2) , x[15] , x[-3]])/1000 #m
        c8b = np.array([8*tlb+x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]/2 , W-(2*tw+x[-3]*3/2) , x[16] , x[-3]])/1000 #m
        c9b = np.array([9*tlb+x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]/2 , W-(2*tw+x[-3]*3/2) , x[17] , x[-3]])/1000 #m

        c1c = np.array([tlc+x[0]/2 , 2*tw+x[-3]*3/2 , x[18] , x[-3]])/1000 #m
        c2c = np.array([2*tlc+x[0]+x[1]/2 , 2*tw+x[-3]*3/2 , x[19] , x[-3]])/1000 #m
        c3c = np.array([3*tlc+x[0]+x[1]+x[2]/2 , 2*tw+x[-3]*3/2 , x[20] , x[-3]])/1000 #m
        c4c = np.array([4*tlc+x[0]+x[1]+x[2]+x[3]/2 ,2*tw+x[-3]*3/2 , x[21] , x[-3]])/1000 #m
        c5c = np.array([5*tlc+x[0]+x[1]+x[2]+x[3]+x[4]/2 , 2*tw+x[-3]*3/2 , x[22] , x[-3]])/1000 #m
        c6c = np.array([6*tlc+x[0]+x[1]+x[2]+x[3]+x[4]+x[5]/2 , 2*tw+x[-3]*3/2 , x[23] , x[-3]])/1000 #m
        c7c = np.array([7*tlc+x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]/2 , 2*tw+x[-3]*3/2 , x[24] , x[-3]])/1000 #m
        c8c = np.array([8*tlc+x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]/2 , 2*tw+x[-3]*3/2 , x[25] , x[-3]])/1000 #m
        c9c = np.array([9*tlc+x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]/2 , 2*tw+x[-3]*3/2 , x[26] , x[-3]])/1000 #m

        c1d = np.array([tld+x[0]/2 , tw+x[-3]/2 , x[27] , x[-3]])/1000 #m
        c2d = np.array([2*tld+x[0]+x[1]/2 , tw+x[-3]/2 , x[28] , x[-3]])/1000 #m
        c3d = np.array([3*tld+x[0]+x[1]+x[2]/2 , tw+x[-3]/2 , x[29] , x[-3]])/1000 #m
        c4d = np.array([4*tld+x[0]+x[1]+x[2]+x[3]/2 , tw+x[-3]/2 , x[30] , x[-3]])/1000 #m
        c5d = np.array([5*tld+x[0]+x[1]+x[2]+x[3]+x[4]/2 , tw+x[-3]/2 , x[31] , x[-3]])/1000 #m
        c6d = np.array([6*tld+x[0]+x[1]+x[2]+x[3]+x[4]+x[5]/2 , tw+x[-3]/2 , x[32] , x[-3]])/1000 #m
        c7d = np.array([7*tld+x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]/2 , tw+x[-3]/2 , x[33] , x[-3]])/1000 #m
        c8d = np.array([8*tld+x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]/2 , tw+x[-3]/2 , x[34] , x[-3]])/1000 #m
        c9d = np.array([9*tld+x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]/2 , tw+x[-3]/2 , x[35] , x[-3]])/1000 #m

        ca = np.concatenate((c1a,c2a,c3a,c4a,c5a,c6a,c7a,c8a,c9a))
        cb = np.concatenate((c1b,c2b,c3b,c4b,c5b,c6b,c7b,c8b,c9b))
        cc = np.concatenate((c1c,c2c,c3c,c4c,c5c,c6c,c7c,c8c,c9c))
        cd = np.concatenate((c1d,c2d,c3d,c4d,c5d,c6d,c7d,c8d,c9d))
        
        # Thermal parameters
        lambda_clay = x[-1]/100 #W/mK

        # Structural parameters
        rho_clay = 2.0841e-6*lambda_clay + 0.55042e-6 #kg/mm3
        resist_clay = 76.766*lambda_clay-26.362 #MPa
        resist_brick = resist_clay*An/Ab #MPa
        resist_brickG2 = 10 #Mpa

        # Calculo U_muro
        if tla < 10 or tlb < 10 or tlc < 10 or tld < 10 or tw < 10:
            U_muro = 50
            if Ah/Ab < 0.5:
                U_muro = 50
            else: U_muro = ThermalAnalysis(lambda_clay,b,ca,cb,cc,cd,dimz,n_cav)
        else:
            U_muro = ThermalAnalysis(lambda_clay,b,ca,cb,cc,cd,dimz,n_cav)

        # Objective 1: minimize weigth
        f1 = rho_clay*H*An #kg

        # Objective 2: minimize wall heat transfer
        f2 = U_muro #W/mK

        # Structural Constraint
        g1 = - resist_brick + resist_brickG2

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
    pop_size = 100, #population size
    n_offsprings = 100, #new individuals from crosssover and mutation (children),
    sampling = IntegerRandomSampling(), #creates initial population
    crossover = SBX(prob=0.9, eta=15, vtype=float, repair=RoundingRepair()),
    mutation = PM(prob=1/n_var, eta=20, vtype=float, repair=RoundingRepair()),
    survival = RankAndCrowdingSurvival(), #Yo Agregue esto
    eliminate_duplicates = True, #element created are different from the ones that already exist
)

#### DEFINE A TERMINATIN CRITERION ####

termination = get_termination("n_gen",30)

#### OPTIMIZE ####

res = minimize(problem,
               algorithm,
               termination,
               seed=1, #ramdom parameter to ensure reproducible results
               save_history=True,
               verbose=True #some printouts during the algorithm´s execution is provided
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



