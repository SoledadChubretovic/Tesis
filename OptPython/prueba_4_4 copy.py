# %% ACTIVAR LICENCIA Y LIBRERIAS
from functions import (
    generate_cavitiy_ansys_parameters
)
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
    BTZ,
    RsiRse,
    λ_glue_mortar,
    prop_ladrillo,
    prop_mort,
    cold_area_index,
    hot_area_index,
    element_size
)
import time
import numpy as np
#from ansys.mapdl.core import launch_mapdl
#mapdl = launch_mapdl()

##### CALCULO TERMICO LADRILLO
start = time.time()
x = [
    # l1, ..., ln
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    # w
    42,
    # W
    118,
    # lambda_clay
    61
]

#### GEOMETRY AND MATERIAL VALUES ####

# vector x is defined as [l1, l2, ..., ln-1, ln, w, W, λ]
_w_pos = len(x) - 3
_W_pos = len(x) - 2
_λ_pos = len(x) - 1

# lower limit fo variables
xl = np.ones(N_VARIABLES - 2) * T_MIN  # mm
xl = np.append(xl, W2_MIN)  # mm
xl = np.append(xl, LAMBDA_CLAY100_MIN)  # mm

# upper limits of variables
xu = np.ones(N_CAVITIES) * l_MAX  # mm
xu = np.append(xu, w_MAX)  # mm
xu = np.append(xu, W2_MAX)  # mm
xu = np.append(xu, LAMBDA_CLAY100_MAX)  # mm

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

print("----------------------------------")

# if there is anything in ANSYS model, finish and clear
mapdl.finish()
mapdl.clear()

# PRE-PROCESSING
mapdl.prep7()

# creates all areas of the brick (perimeter and cavities)
brick_perimeter = mapdl.blc5(*b)
for row in matrix:
    mapdl.blc5(*row)

# substract all other areas from brick_perimeter
# creates a new area with cavities
# the index of this area is the number of cavities + 2
# erases all other areas
cavities = mapdl.asba(brick_perimeter, 'all')

# plot final area with its index
mapdl.allsel()
#mapdl.aplot(show_area_numbering=True)

# extrude area to create a volume
mapdl.vext(cavities, dz = DIM_Z)  # m

# plot volume showing index numbers of each area
#mapdl.vplot(show_area_numbering=True)

# material properties 
# thermal conductivity for isotropic material in W/mK
Id_arcilla = 1
mapdl.mp("KXX", Id_arcilla, λ_clay) # W/mK
mapdl.mp("KYY", Id_arcilla, λ_clay) # W/mK
mapdl.mp("KZZ", Id_arcilla, λ_clay) # W/mK

# finite element
# assign element type for the analysis
# SOLID70: 3D thermal element with 8 nodes 
# 1 degree of freedom per node: temperature
Id_ladrillo = 1
mapdl.et(Id_ladrillo, "Solid70")

# asigne attributes to the volume
# real constat = 0
mapdl.vatt(Id_arcilla, 0, Id_ladrillo)

# mesh
# mesh shape
mapdl.mshape(1, "3D") # 1: tetrahedral shape
# use main tetrahedral mesher
mapdl.mopt("vmesh", "default")
# size of finite elements
mapdl.esize(element_size)
# permit that elements expand when advancing to the interior of the volume
# guide the mesh froma  fine mesh on the boundry to a coarse mesh on the exterior
#mapdl.mopt("expnd", 3)
# controlling mesh transitioning
#mapdl.mopt("trans", 1.5) #
# mesh all volumes (in this case we have only one)
mapdl.vmesh("all")
#mapdl.eplot() # plot meshed volume

# border conditions and DOF constraints
# select desired areas and asign temperatures
# temperatures assigned to creat heat flux
Thot = 30  # C
Tcold = 18  # C
mapdl.asel("S", vmin = cold_area_index)
mapdl.da("all", "temp", Tcold)
mapdl.asel("S", vmin = hot_area_index)
mapdl.da("all", "temp", Thot)
out = mapdl.allsel()

# SOLVE
# solve model
mapdl.run("/SOLU")
mapdl.solve()
out = mapdl.finish()

# POST-PROCESSING
# obtain results
mapdl.post1()
mapdl.set("last", "last")

# plot brick with temperatures
# mapdl.post_processing.plot_nodal_temperature()

# calculate heat flux and U_value (wall thermal transmitance)

# cold area (area with 18C degrees forced)
mapdl.set("last", "last")
mapdl.asel("S", vmin = cold_area_index)
nodes = mapdl.nsla()
# first node of cold area (minimum index)
min_nodenum_cold = int(mapdl.get("min_nodenum_cold", "node", "0", "num", "min"))
# last node of cold area (maximum index)
max_nodenum_cold = int(mapdl.get("max_nodenum_cold", "node", "0", "num", "max"))
# number of selected nodes (all nodes of cold area)
nb_selected_nodes_cold = mapdl.mesh.n_node 

# fcold is a vector that stores the thermal flux in y axis for all nodes in cold area
j = 0
fcold = np.zeros(nb_selected_nodes_cold)
for i in range(min_nodenum_cold, max_nodenum_cold + 1):
    fcold[j] = mapdl.get("fcold", "node", i, "tf", "y")
    j = j + 1
# mean flux in cold area in W/m2
fcold = sum(abs(fcold))/len(fcold)

# hot area (area with 30C degrees forced)
mapdl.set("last", "last")
mapdl.asel("S", vmin = hot_area_index)
nodes = mapdl.nsla()
# first node of hot area (minimum index)
min_nodenum_hot = int(mapdl.get("min_nodenum_hot", "node", "0", "num", "min"))
# las node of hot area (maximum index)
max_nodenum_hot = int(mapdl.get("max_nodenum_hot", "node", "0", "num", "max"))
# number of selected nodes (all nodes of hot area)
nb_selected_nodes_hot = mapdl.mesh.n_node  # numero de nodos seleccionados

# fhot is a vector that stores the thermal flux in y axis for all nodes in hot area
j = 0
fhot = np.zeros(nb_selected_nodes_hot)
for i in range(min_nodenum_hot, max_nodenum_hot + 1):
    fhot[j] = mapdl.get("fhot", "node", i, "tf", "y")
    j = j + 1
# mean flux in hot area in W/m2
fhot = sum(abs(fhot)) / len(fhot)

# heat flux per area in y axis in W/m2
flux_area = (fcold + fhot)/2  # W/m2

# equivalent thermal conductivity in W/mK
λeq = flux_area * b[3] / (Thot - Tcold)  # W/mK

# wall thermal conductance in W/m2K
C_muro = (prop_ladrillo * λeq + prop_mort * λ_glue_mortar) / b[3]  # W/m2K
# wall thermal transmitance in W/m2K
U_muro = 1 / (RsiRse + 1 / C_muro)  # W/m2K

print(U_muro)

#mapdl.post_processing.plot_nodal_temperature()

# convection iterations
# until here the program considers only convection
# convection and radiation in cavities is consider in H paramter from ISO6946
# iterations are needed because H parameter depends of cavity temperature
# wich also varies with H value (recursive)
# 3 iterations are enough to converge in a valid result
iteracion = 0
while iteracion < 3:
    
    print(iteracion)
    # calculus of mean temperature of each cavity
    # areas_convection: list with all index for inner areas
    areas_convection = []
    for i in range(6, N_CAVITIES + 2):
        areas_convection = np.append(areas_convection, i)
    for i in range(N_CAVITIES + 3, N_CAVITIES * 4 + 7):
        areas_convection = np.append(areas_convection, i)

    T_conv = np.zeros(len(areas_convection))
    for idx, i in enumerate(areas_convection):
        mapdl.asel("S", vmin=i)
        mapdl.nsla()
        # T: temperatures for each node for each inner area
        T = mapdl.post_processing.nodal_temperature()
        T_conv[idx] = np.mean(T)  # mean temperature for each inner area

    Tm_conv = np.zeros(N_CAVITIES)  # mean temperature for each cavity
    for idx, i in enumerate(range(0, len(T_conv), 4)):
        Tm_conv[idx] = np.mean(T_conv[i:i+4])

    # H parameter calculus (W/m2K)
    # this parameter considers convection and radiation for small cavities
    # based ISO6946

    # d: dimension of cavity in direction of the heat flux
    # d = w
    d = np.zeros(N_CAVITIES)
    for i in range(0, len(d)):
        d[i] = matrix[i, 3]  # m
    
    # b1: dimension of cavity perpendicular to the heat flux
    # b1 = l
    b1 = np.zeros(N_CAVITIES)
    for i in range(0, len(b1)):
        b1[i] = matrix[i, 2]  # m

    # hr0 = 4*btz*Tmn^3
    hr0 = 4 * BTZ * (np.array(Tm_conv + 273.15) ** 3)

    # ha value varies with value of d
    # ha = max(0.025 / d , 1.25)
    ha = np.zeros(N_CAVITIES)
    for i in range(0, len(ha)):
        ha[i] = max(0.025 / d[i], 1.25)

    # hr: huge formula
    emissivities = (1 / RAD) + (1 / RAD)
    square_root = np.sqrt(1 + (d ** 2) / (b1 ** 2))
    last_part = 2 / (1 + square_root - d / b1)
    hr = hr0 / (emissivities - 2 + last_part)

    # convection and radiation parameter
    # H_conv is a vector with the H parameter for each cavity
    H_conv = hr + ha  # W/m2K

    # surface load convection and radiation
    mapdl.run("/SOLU")
    mapdl.sfcum("all", "REPL")  # replace load in area

    # asign H value in each 4 areas of each cavity
    j = 0
    for i in range(0, N_CAVITIES):
        # cuadrupleta stores the first index and last index of the 4 areas of a cavity
        cuadrupleta = [i*4, i*4+3]
        # select areas of a cavity
        mapdl.asel("S", vmin = areas_convection[cuadrupleta[0]], vmax = areas_convection[cuadrupleta[-1]])
        mapdl.nsla()  # select nodes in selected areas
        mapdl.sf("all", "conv", H_conv[j], Tm_conv[j])
        j = j + 1
    out = mapdl.allsel()

    # SOLVE
    mapdl.solve()
    out = mapdl.finish()

    # POST-PROCESSING
    mapdl.post1()
    mapdl.set("last", "last")

    # plot brick with temperatures
    # mapdl.post_processing.plot_nodal_temperature()

    iteracion = iteracion + 1

# plot brick with temperatures
# mapdl.post_processing.plot_nodal_temperature()

# calculate heat flux and U_value (wall thermal transmitance)

# cold area (area with 18C degrees forced)
mapdl.set("last", "last")
mapdl.asel("S", vmin = cold_area_index)
nodes = mapdl.nsla()
# first node of cold area (minimum index)
min_nodenum_cold = int(mapdl.get("min_nodenum_cold", "node", "0", "num", "min"))
# last node of cold area (maximum index)
max_nodenum_cold = int(mapdl.get("max_nodenum_cold", "node", "0", "num", "max"))
# number of selected nodes (all nodes of cold area)
nb_selected_nodes_cold = mapdl.mesh.n_node 

# fcold is a vector that stores the thermal flux in y axis for all nodes in cold area
j = 0
fcold = np.zeros(nb_selected_nodes_cold)
for i in range(min_nodenum_cold, max_nodenum_cold + 1):
    fcold[j] = mapdl.get("fcold", "node", i, "tf", "y")
    j = j + 1
# mean flux in cold area in W/m2
fcold = sum(abs(fcold))/len(fcold)

# hot area (area with 30C degrees forced)
mapdl.set("last", "last")
mapdl.asel("S", vmin = hot_area_index)
nodes = mapdl.nsla()
# first node of hot area (minimum index)
min_nodenum_hot = int(mapdl.get("min_nodenum_hot", "node", "0", "num", "min"))
# las node of hot area (maximum index)
max_nodenum_hot = int(mapdl.get("max_nodenum_hot", "node", "0", "num", "max"))
# number of selected nodes (all nodes of hot area)
nb_selected_nodes_hot = mapdl.mesh.n_node  # numero de nodos seleccionados

# fhot is a vector that stores the thermal flux in y axis for all nodes in hot area
j = 0
fhot = np.zeros(nb_selected_nodes_hot)
for i in range(min_nodenum_hot, max_nodenum_hot + 1):
    fhot[j] = mapdl.get("fhot", "node", i, "tf", "y")
    j = j + 1
# mean flux in hot area in W/m2
fhot = sum(abs(fhot)) / len(fhot)

# heat flux per area in y axis in W/m2
flux_area = (fcold + fhot)/2  # W/m2

# equivalent thermal conductivity in W/mK
λeq = flux_area * b[3] / (Thot - Tcold)  # W/mK

# wall thermal conductance in W/m2K
C_muro = (prop_ladrillo * λeq + prop_mort * λ_glue_mortar) / b[3]  # W/m2K
# wall thermal transmitance in W/m2K
U_muro = 1 / (RsiRse + 1 / C_muro)  # W/m2K

print(U_muro)

end = time.time()
print("Time: ", (end-start), "s")

# %%
