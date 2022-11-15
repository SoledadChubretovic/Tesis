# %% ACTIVAR LICENCIA Y LIBRERIAS
from functions import (
    generate_cavitiy_ansys_parameters
)

from k_constants import (
    N_CAVITIES,
    DIM_Z,
    RAD,
    BTZ,
    RsiRse,
    λ_GLUE_MORTAR,
    PROP_BRICK,
    prop_mortar,
    COLD_AREA_INDEX,
    HOT_AREA_INDEX,
    ELEMENT_SIZE,
    T_HOT,
    T_COLD,
    U_MURO_INVALID
)

import numpy as np

# start licence instance 
# this is already done in OpenANSYS
# from ansys.mapdl.core import launch_mapdl
# mapdl = launch_mapdl()

def ThermalAnalysis(b, matrix, λ_clay):
    try:
        # thermal analysis of the brick
        # vector x is defined as [l1, l2, ..., ln-1, ln, w, W, λ]

        # if there is anything in ANSYS model, finish and clear
        mapdl.finish()
        mapdl.clear()

        # PRE-PROCESSING
        mapdl.prep7()

        # define title for analysis
        mapdl.title("Tesis")

        # creates all areas of the brick (perimeter and cavities)
        brick_perimeter = mapdl.blc5(*b)
        for row in matrix:
            mapdl.blc5(*row)

        # substract all other areas from brick_perimeter
        # creates a new area with cavities
        # the index of this area is the number of cavities + 2
        # erases all other areas
        cavities = mapdl.asba(brick_perimeter, 'all')

        # plot final area with indices
        # mapdl.allsel()
        # mapdl.aplot(show_area_numbering=True)

        # extrude area to create a volume
        mapdl.allsel()
        mapdl.vext(cavities, dz = DIM_Z)  # m

        # plot volume showing index numbers of each area
        # mapdl.vplot(show_area_numbering = True)

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

        # real constant
        # real constat set to 0 because there is not mass movement
        Id_real_constant = 1
        mapdl.r(Id_real_constant,0)

        # asigne attributes to the volume
        mapdl.vatt(Id_arcilla, Id_real_constant, Id_ladrillo)

        # mesh
        # mesh shape
        mapdl.mshape(1, "3D") # 1: tetrahedral shape
        # size of finite elements
        mapdl.esize(ELEMENT_SIZE)
        # permit that elements expand when advancing to the interior of the volume
        # guide the mesh from a fine mesh on the boundry to a coarse mesh on the exterior
        # mapdl.mopt("expnd", 1.2)
        # controlling mesh transitioning
        # mapdl.mopt("trans", 1.2) #
        # use main tetrahedral mesher
        mapdl.mopt("vmesh", "default")
        # mesh all volumes (in this case we have only one)
        mapdl.vmesh("all")
        # plot meshed volume
        # mapdl.eplot()

        # border conditions and DOF constraints
        # select desired areas and asign temperatures
        # temperatures assigned to creat heat flux
        Thot = T_HOT  # C
        Tcold = T_COLD  # C
        mapdl.asel("S", vmin = COLD_AREA_INDEX)
        mapdl.da("all", "temp", Tcold)
        mapdl.asel("S", vmin = HOT_AREA_INDEX)
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

        # plot brick with temperatures
        # mapdl.post_processing.plot_nodal_temperature()

        # convection iterations
        # until here the program considers only conduction
        # convection and radiation in cavities are considered in H parameter from ISO6946
        # iterations are needed because H parameter depends of cavity temperature
        # which also varies with H value
        # in this case 3 iterations are enough to converge in a valid result

        iteracion = 0
        while iteracion < 3:

            # calculation of mean temperature of each cavity
            # areas_convection: list with all index for inner areas
            areas_convection = []
            for i in range(6, N_CAVITIES + 2):
                areas_convection = np.append(areas_convection, i)
            for i in range(N_CAVITIES + 3, N_CAVITIES * 4 + 7):
                areas_convection = np.append(areas_convection, i)

            # this loop gets the temperature of each node of an area 
            # and calculates the mean value
            # T_conv: mean temperature for each cavity
            T_conv = np.zeros(len(areas_convection))
            for idx, i in enumerate(areas_convection):
                mapdl.asel("S", vmin = i)
                mapdl.nsla()
                # T: temperatures for each node for each inner area
                T = mapdl.post_processing.nodal_temperature()
                T_conv[idx] = np.mean(T)

            # this loop calculates the mean temperature of each cavity 
            # with temperatures of each area of the cavity
            # Tm_conv: mean temperature for each cavity
            Tm_conv = np.zeros(N_CAVITIES)
            for idx, i in enumerate(range(0, len(T_conv), 4)):
                Tm_conv[idx] = np.mean(T_conv[i:i+4])

            # H parameter calculation (W/m2K)
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

            # hr: huge formula for radiation parameter
            emissivities = (1 / RAD) + (1 / RAD)
            square_root = np.sqrt(1 + (d ** 2) / (b1 ** 2))
            last_part = 2 / (1 + square_root - d / b1)
            hr = hr0 / (emissivities - 2 + last_part)

            # convection and radiation parameter
            # H_conv is a vector with the H parameter for each cavity
            H_conv = hr + ha  # W/m2K

            # surface load convection and radiation
            mapdl.run("/SOLU")
            # replace loads in areas
            mapdl.sfcum("all", "REPL")

            # assign H value in each 4 areas of each cavity
            for i in range(0, N_CAVITIES):
                cuadrupleta = [i*4, i*4+3]
                # select areas
                mapdl.asel("S", "area", vmin = areas_convection[cuadrupleta[0]], vmax = areas_convection[cuadrupleta[-1]])
                mapdl.nsla("S", 1)  # select nodes from selected areas
                mapdl.sf("all", "conv", H_conv[i], Tm_conv[i]) #assign H and Tm

            out = mapdl.allsel()

            # SOLVE
            mapdl.solve()
            out = mapdl.finish()

            # POST-PROCESSING
            mapdl.post1()
            # obtain results from the last step
            mapdl.set("last", "last")

            # plot brick with temperatures
            # mapdl.post_processing.plot_nodal_temperature()

            iteracion = iteracion + 1

        # calculate heat flux and U_value (wall thermal transmitance)

        # cold area (area with 18C degrees forced)
        mapdl.set("last", "last")
        # select cold area and all its nodes
        mapdl.asel("S", vmin = COLD_AREA_INDEX)
        mapdl.nsla()
        # first node of cold area (minimum index)
        min_nodenum_cold = int(mapdl.get("min_nodenum_cold", "node", "0", "num", "min"))
        # last node of cold area (maximum index)
        max_nodenum_cold = int(mapdl.get("max_nodenum_cold", "node", "0", "num", "max"))
        # number of selected nodes (all nodes of cold area)
        nb_selected_nodes_cold = mapdl.mesh.n_node 

        # fcold: vector with thermal flux in y axis for all nodes in cold area
        j = 0
        fcold = np.zeros(nb_selected_nodes_cold)
        for i in range(min_nodenum_cold, max_nodenum_cold + 1):
            fcold[j] = mapdl.get("fcold", "node", i, "tf", "y")
            j = j + 1
        # mean flux in cold area in W/m2
        fcold = sum(abs(fcold))/len(fcold)

        # hot area (area with 30C degrees forced)
        mapdl.set("last", "last")
        # select cold area and all its nodes
        mapdl.asel("S", vmin = HOT_AREA_INDEX)
        mapdl.nsla()
        # first node of hot area (minimum index)
        min_nodenum_hot = int(mapdl.get("min_nodenum_hot", "node", "0", "num", "min"))
        # las node of hot area (maximum index)
        max_nodenum_hot = int(mapdl.get("max_nodenum_hot", "node", "0", "num", "max"))
        # number of selected nodes (all nodes of hot area)
        nb_selected_nodes_hot = mapdl.mesh.n_node

        # fhot: vector with thermal flux in y axis for all nodes in hot area
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
        C_muro = (PROP_BRICK * λeq + prop_mortar * λ_GLUE_MORTAR) / b[3]  # W/m2K
        # wall thermal transmitance in W/m2K
        U_muro = 1 / (RsiRse + 1 / C_muro)  # W/m2K

        f = open(r"C:\Users\nchubretovic\OneDrive - Entel\Escritorio\Sole\Tesis\OptPython\1_Model\Error.txt", "a")
        f.write("Umuro = " + U_muro + "W/m2K" + "\n")
        f.close
    except:
        U_muro = U_MURO_INVALID
        f.write("Fail --- " + " Error: " + "\n")
        for i in range(len(sys.exc_info())):
            f.write(str(sys.exc_info()[i]) + "\n")
        f.colse()

    return U_muro
# %%
