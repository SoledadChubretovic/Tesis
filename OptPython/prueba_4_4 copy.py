#%% ACTIVAR LICENCIA Y LIBRERIAS
from ansys.mapdl.core import launch_mapdl
mapdl = launch_mapdl()

# %% CALCULO TERMICO LADRILLO

import numpy as np
import time

from k_constants import (
    CAVITIES_PER_ROW,
    N_CAVITIES,
    N_VARIABLES,
    L,
    H,
    DIM_Z,
    T_MIN,
    W2_MIN,
    W2_MAX,
    LAMBDA_CLAY100_MIN,
    LAMBDA_CLAY100_MAX,
    l_MAX,
    w_MAX,
    RAD,
    BTZ
)
from RH import (
    generate_cavitiy_ansys_parameters
)

x = [
    # l1, ..., ln
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    # w
    42,
    # W
    118,
    # lambda_clay
    61
]

# lower limit fo variables
xl = np.ones(N_VARIABLES - 2) * T_MIN # mm
xl = np.append(xl, W2_MIN) # mm
xl = np.append(xl, LAMBDA_CLAY100_MIN) # mm

# upper limits of variables
xu = np.ones(N_CAVITIES) * l_MAX # mm
xu = np.append(xu, w_MAX) # mm
xu = np.append(xu, W2_MAX) # mm
xu = np.append(xu, LAMBDA_CLAY100_MAX) # mm

#### GEOMETRY AND MATERIAL VALUES ####

# vector is defined as [l1, l2, ..., ln-1, ln, w, W, λ]
_w_pos = len(x) - 3
_W_pos = len(x) - 2
_λ_pos = len(x) - 1

# width of the brick
W = x[_W_pos] * 2 # mm

# width of every cavity
w = x[_w_pos] #mm

# area of the cavities within the brick ()
Ah = x[_w_pos] * sum(x[:N_CAVITIES])
# total area of the brick
Ab = W * L
# net area of the brick (total - cavities)
An = Ab - Ah

# [centro x, centro y, largo x, largo y]
b = np.array([L/2,W/2,L,W])/1000 # m
matrix = generate_cavitiy_ansys_parameters(x, W, w) / 1000 # m

# Thermal parameters
λ_clay = x[_λ_pos] / 100 # W/mK

elemento = [0.006]

for elefini in elemento:
    print("----------------------------------")
    print(elefini)
    start = time.time()

    mapdl.finish()
    mapdl.clear()
    mapdl.prep7()

    ###### GEOMETRY ######
    brick_perimeter = mapdl.blc5(*b)
    for row in matrix:
        mapdl.blc5(*row)

    cavities = mapdl.asba(brick_perimeter, 'all')
    mapdl.allsel()
    mapdl.aplot(show_area_numbering=True)
    mapdl.vext(cavities, dz = DIM_Z) # m
    mapdl.vplot(show_area_numbering=True)

    ###### ASIGNACION DE ATRIBUTOS ###### 

    # Material (Conductividad termica, material isotrópico)
    Id_arcilla = 1
    mapdl.mp("KXX", Id_arcilla, λ_clay)
    mapdl.mp("KYY", Id_arcilla, λ_clay)
    mapdl.mp("KZZ", Id_arcilla, λ_clay)

    # Elemento finito
    Id_ladrillo = 1
    mapdl.et(Id_ladrillo, "Solid70") # element reference p.297 "SOLID70"

    # Asignar todo a ladrillo
    mapdl.vatt(Id_arcilla, 0, Id_ladrillo)

    ###### MALLADO ######
    mapdl.mshape(1, "3D")
    mapdl.esize(elefini)
    mapdl.mopt("vmesh", "default")
    mapdl.mopt("expnd", 3)
    mapdl.mopt("trans", 1.5)
    mapdl.vmesh("all")
    mapdl.eplot()

    ###### CONDICIONES DE BORDE DOF CONSTRAINS ######

    Thot = 30 # C
    Tcold = 18 # C

    mapdl.asel("S", vmin=2)
    mapdl.da("all", "temp", Tcold)
    mapdl.asel("S", vmin=4)
    mapdl.da("all", "temp", Thot)
    out = mapdl.allsel()

    ###### SOLVE ###### 
    mapdl.run("/SOLU")
    mapdl.solve()
    out = mapdl.finish()

    # POST-PROCESSING
    mapdl.post1()
    mapdl.set("last","last")
    mapdl.post_processing.plot_nodal_temperature()

    ###### FLUJO DE CALOR ###### 

    # Cara fria
    mapdl.set("last","last")
    mapdl.asel("S", vmin=2)
    nodes = mapdl.nsla()

    min_nodenum_cold = int(mapdl.get("min_nodenum_cold","node","0","num","min")) #numero minimo de nodo seleccionado
    max_nodenum_cold = int(mapdl.get("max_nodenum_cold","node","0","num","max")) #numero maximo de nodo seleccionado
    nb_selected_nodes_cold = mapdl.mesh.n_node #numero de nodos seleccionados

    j = 0
    fcold = np.zeros(nb_selected_nodes_cold) 
    for i in range(min_nodenum_cold,max_nodenum_cold + 1):
        fcold[j] = mapdl.get("fcold","node",i,"tf","y")
        j = j + 1

    fcold = sum(abs(fcold))/len(fcold)

    # Cara caliente
    mapdl.set("last","last")
    mapdl.asel("S", vmin=4)
    nodes = mapdl.nsla()

    min_nodenum_hot = int(mapdl.get("min_nodenum_hot","node","0","num","min")) #numero minimo de nodo seleccionado
    max_nodenum_hot = int(mapdl.get("max_nodenum_hot","node","0","num","max")) #numero maximo de nodo seleccionado
    nb_selected_nodes_hot = mapdl.mesh.n_node #numero de nodos seleccionados

    j = 0
    fhot = np.zeros(nb_selected_nodes_hot) 
    for i in range(min_nodenum_hot,max_nodenum_hot + 1):
        fhot[j] = mapdl.get("fhot","node",i,"tf","y")
        j = j + 1

    fhot = sum(abs(fhot))/len(fhot)
    flux_area = (fcold+fhot)/2 #W/m2

    ###### CONDUCTIVIDAD EQUIVALENTE ###### 

    keq = flux_area*b[3]/(Thot-Tcold) #W/mK

    ###### TRANSMITANCIA DEL MURO ######

    RsiRse = 0.17 #m2K/W
    k_mort_pega = 0.23 #W/mK
    prop_ladrillo = 0.84
    prop_mort = 1 - prop_ladrillo
    C_muro = (prop_ladrillo*keq+prop_mort*k_mort_pega)/b[3] #W/m2K
    U_muro = 1/(RsiRse+1/C_muro) #W/m2K

    print(U_muro)

    ###### ITERACIONES CONVECCION ###### 

    iteracion = 0

    while iteracion < 3:
        print(iteracion)
        # Calculo Tm
        areas_convection = [] # lista con todas las areas interiores (numeros correctos, saltandose la que corresponde)
        for i in range(6, N_CAVITIES + 2):
            areas_convection = np.append(areas_convection, i)
        for i in range(N_CAVITIES + 3, N_CAVITIES * 4 + 7):
            areas_convection = np.append(areas_convection, i)

        T_conv = np.zeros(len(areas_convection))
        for idx, i in enumerate(areas_convection):
            mapdl.asel("S", vmin=i)
            mapdl.nsla()
            T = mapdl.post_processing.nodal_temperature() #temperaturas en cada nodo del area
            T_conv[idx]= np.mean(T) #temperatura promedio del area

        Tm_conv = np.zeros(N_CAVITIES) #temperatura media para cada cavidad
        for idx, i in enumerate(range(0,len(T_conv),4)):
            Tm_conv[idx] = np.mean(T_conv[i:i+4])

        # Calculo parametro H (W/m2K)

        d = np.zeros(N_CAVITIES)
        for i in range(0,len(d)):
            d[i]= matriz[i,3] #m  

        b1 = np.zeros(N_CAVITIES)
        for i in range(0,len(b1)):
            b1[i]= matriz[i,2] #m

        hr0 = 4*BTZ*(np.array(Tm_conv+273.15)**3)

        ha = np.zeros(N_CAVITIES)
        for i in range(0,len(ha)):
            ha[i] = max(0.025/d[i],1.25)

        hr = hr0/((1/RAD)+(1/RAD)-2+(2/(1+np.sqrt(1+(d**2)/(b1**2))-d/b1)))

        H_conv = hr + ha #W/m2K

        # SURFACE LOAD CONVECTION

        mapdl.run("/SOLU")
        mapdl.sfcum("all","REPL") #reemplaza las cargas en areas

        j = 0
        for i in range(0,int(len(a_conv)/4)):
            cuadrupleta = [i*4,i*4+3]
            mapdl.asel("S", vmin=a_conv[cuadrupleta[0]],vmax=a_conv[cuadrupleta[-1]]) #selecciona areas
            mapdl.nsla() #selecciona nodos asociados a las areas seleccionadas
            mapdl.sf("all","conv",H_conv[j],Tm_conv[j])
            j = j+1

        out = mapdl.allsel()

        # SOLVE
        mapdl.solve()
        out = mapdl.finish()

        # POST-PROCESSING
        mapdl.post1()
        mapdl.set("last","last")
        #mapdl.post_processing.plot_nodal_temperature()
        
        iteracion = iteracion + 1

        ###### FLUJO DE CALOR ###### 

        # Cara fria
        mapdl.set("last","last")
        mapdl.asel("S", vmin=2)
        nodes = mapdl.nsla()

        min_nodenum_cold = int(mapdl.get("min_nodenum_cold","node","0","num","min")) #numero minimo de nodo seleccionado
        max_nodenum_cold = int(mapdl.get("max_nodenum_cold","node","0","num","max")) #numero maximo de nodo seleccionado
        nb_selected_nodes_cold = mapdl.mesh.n_node #numero de nodos seleccionados

        j = 0
        fcold = np.zeros(nb_selected_nodes_cold) 
        for i in range(min_nodenum_cold,max_nodenum_cold + 1):
            fcold[j] = mapdl.get("fcold","node",i,"tf","y")
            j = j + 1

        fcold = sum(abs(fcold))/len(fcold)

        # Cara caliente
        mapdl.set("last","last")
        mapdl.asel("S", vmin=4)
        nodes = mapdl.nsla()

        min_nodenum_hot = int(mapdl.get("min_nodenum_hot","node","0","num","min")) #numero minimo de nodo seleccionado
        max_nodenum_hot = int(mapdl.get("max_nodenum_hot","node","0","num","max")) #numero maximo de nodo seleccionado
        nb_selected_nodes_hot = mapdl.mesh.n_node #numero de nodos seleccionados

        j = 0
        fhot = np.zeros(nb_selected_nodes_hot) 
        for i in range(min_nodenum_hot,max_nodenum_hot + 1):
            fhot[j] = mapdl.get("fhot","node",i,"tf","y")
            j = j + 1

        fhot = sum(abs(fhot))/len(fhot)
        flux_area = (fcold+fhot)/2 #W/m2

        ###### CONDUCTIVIDAD EQUIVALENTE ###### 

        keq = flux_area*b[3]/(Thot-Tcold) #W/mK

        ###### TRANSMITANCIA DEL MURO ######

        RsiRse = 0.17 #m2K/W
        k_mort_pega = 0.23 #W/mK
        prop_ladrillo = 0.84
        prop_mort = 1 - prop_ladrillo
        C_muro = (prop_ladrillo*keq+prop_mort*k_mort_pega)/b[3] #W/m2K
        U_muro = 1/(RsiRse+1/C_muro) #W/m2K

        print(U_muro)

    end = time.time()
    print("Time: ", (end-start), "s")
#%%otro

# ###### FLUJO DE CALOR ###### 

# # Cara fria
# mapdl.set("last","last")
# mapdl.asel("S", vmin=2)
# nodes = mapdl.nsla()

# min_nodenum_cold = int(mapdl.get("min_nodenum_cold","node","0","num","min")) #numero minimo de nodo seleccionado
# max_nodenum_cold = int(mapdl.get("max_nodenum_cold","node","0","num","max")) #numero maximo de nodo seleccionado
# nb_selected_nodes_cold = mapdl.mesh.n_node #numero de nodos seleccionados

# j = 0
# fcold = np.zeros(nb_selected_nodes_cold) 
# for i in range(min_nodenum_cold,max_nodenum_cold + 1):
#     fcold[j] = mapdl.get("fcold","node",i,"tf","y")
#     j = j + 1

# fcold = sum(abs(fcold))/len(fcold)

# # Cara caliente
# mapdl.set("last","last")
# mapdl.asel("S", vmin=4)
# nodes = mapdl.nsla()

# min_nodenum_hot = int(mapdl.get("min_nodenum_hot","node","0","num","min")) #numero minimo de nodo seleccionado
# max_nodenum_hot = int(mapdl.get("max_nodenum_hot","node","0","num","max")) #numero maximo de nodo seleccionado
# nb_selected_nodes_hot = mapdl.mesh.n_node #numero de nodos seleccionados

# j = 0
# fhot = np.zeros(nb_selected_nodes_hot) 
# for i in range(min_nodenum_hot,max_nodenum_hot + 1):
#     fhot[j] = mapdl.get("fhot","node",i,"tf","y")
#     j = j + 1

# fhot = sum(abs(fhot))/len(fhot)
# flux_area = (fcold+fhot)/2 #W/m2

# ###### CONDUCTIVIDAD EQUIVALENTE ###### 

# keq = flux_area*b[3]/(Thot-Tcold) #W/mK

# ###### TRANSMITANCIA DEL MURO ######

# RsiRse = 0.17 #m2K/W
# k_mort_pega = 0.23 #W/mK
# prop_ladrillo = 0.84
# prop_mort = 1 - prop_ladrillo
# C_muro = (prop_ladrillo*keq+prop_mort*k_mort_pega)/b[3] #W/m2K
# U_muro = 1/(RsiRse+1/C_muro) #W/m2K

# print(U_muro)
#end = time.time()
#print("Time: ", (end-start), "s")

# %%
