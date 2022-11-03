#%% ACTIVAR LICENCIA Y LIBRERIAS

# Unidades: m,W,C
#from ansys.mapdl.core import launch_mapdl
#mapdl = launch_mapdl()

# %% CALCULO TERMICO LADRILLO
import numpy as np
import time

 #### GEOMETRY AND MATERIAL VALUES ####
# Number of variables, constrains and objective functions
n_cav_L = 4
n_cav_W = 4
n_cav = n_cav_L*n_cav_W
n_var = n_cav + 3 #number of variables (li,w,W,lambda_clay)

# Upper and lower limit
L = 154 #mm
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


x = [15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15, 42, 118, 61]

#### GEOMETRY AND MATERIAL VALUES ####
H = 113 #mm
W = x[-2]*2 #mm

# Areas
Ah = x[-3]*(sum(x)-x[-1]-x[-2]-x[-3])
Ab = W*L
An = Ab-Ah

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

count = 0
c1a = np.array([tla+x[count]/2 , W-(tw+x[-3]/2) , x[count] , x[-3]])/1000 #m
c2a = np.array([2*tla+x[count]+x[count+1]/2 , W-(tw+x[-3]/2) , x[count+1] , x[-3]])/1000 #m
c3a = np.array([3*tla+x[count]+x[count+1]+x[count+2]/2 , W-(tw+x[-3]/2) , x[count+2] , x[-3]])/1000 #m
c4a = np.array([4*tla+x[count]+x[count+1]+x[count+2]+x[count+3]/2 , W-(tw+x[-3]/2) , x[count+3] , x[-3]])/1000 #m

count = count + n_cav_L
c1b = np.array([tlb+x[count]/2 , W-(2*tw+x[-3]*3/2) , x[count] , x[-3]])/1000 #m
c2b = np.array([2*tlb+x[count]+x[count+1]/2 , W-(2*tw+x[-3]*3/2) , x[count+1] , x[-3]])/1000 #m
c3b = np.array([3*tlb+x[count]+x[count+1]+x[count+2]/2 , W-(2*tw+x[-3]*3/2) , x[count+2] , x[-3]])/1000 #m
c4b = np.array([4*tlb+x[count]+x[count+1]+x[count+2]+x[count+3]/2 , W-(2*tw+x[-3]*3/2) , x[count+3] , x[-3]])/1000 #m

count = count + n_cav_L
c1c = np.array([tlc+x[count]/2 , 2*tw+x[-3]*3/2 , x[count] , x[-3]])/1000 #m
c2c = np.array([2*tlc+x[count]+x[count+1]/2 , 2*tw+x[-3]*3/2 , x[count+1] , x[-3]])/1000 #m
c3c = np.array([3*tlc+x[count]+x[count+1]+x[count+2]/2 , 2*tw+x[-3]*3/2 , x[count+2] , x[-3]])/1000 #m
c4c = np.array([4*tlc+x[count]+x[count+1]+x[count+2]+x[count+3]/2 ,2*tw+x[-3]*3/2 , x[count+3] , x[-3]])/1000 #m

count = count + n_cav_L
c1d = np.array([tld+x[count]/2 , tw+x[-3]/2 , x[count] , x[-3]])/1000 #m
c2d = np.array([2*tld+x[count]+x[count+1]/2 , tw+x[-3]/2 , x[count+1] , x[-3]])/1000 #m
c3d = np.array([3*tld+x[count]+x[count+1]+x[count+2]/2 , tw+x[-3]/2 , x[count+2] , x[-3]])/1000 #m
c4d = np.array([4*tld+x[count]+x[count+1]+x[count+2]+x[count+3]/2 , tw+x[-3]/2 , x[count+3] , x[-3]])/1000 #m


ca = np.concatenate((c1a,c2a,c3a,c4a))
cb = np.concatenate((c1b,c2b,c3b,c4b))
cc = np.concatenate((c1c,c2c,c3c,c4c))
cd = np.concatenate((c1d,c2d,c3d,c4d))

# Thermal parameters
lambda_clay = x[-1]/100 #W/mK

matriz = np.array([c1a,c2a,c3a,c4a,c1b,c2b,c3b,c4b,c1c,c2c,c3c,c4c,c1d,c2d,c3d,c4d])

#elemento = [0.0013,0.0012,0.0011,0.00115,0.00112,0.00111]
elemento = [0.003]

for elefini in elemento:
    print("----------------------------------")
    print(elefini)
    start = time.time()

    mapdl.finish()
    mapdl.clear()
    mapdl.prep7()

    ###### GEOMETRIA ######
    brick_perimeter = mapdl.blc5(b[0],b[1],b[2],b[3])
    for row in matriz:
        mapdl.blc5(row[0],row[1],row[2],row[3])

    cavities = mapdl.asba(brick_perimeter, 'all')
    mapdl.allsel()
    #mapdl.aplot(show_area_numbering=True)
    mapdl.vext(cavities,dz = dimz) #m
    #mapdl.aplot(show_area_numbering=True)


    ###### ASIGNACION DE ATRIBUTOS ###### 

    # Material (Conductividad termica, material isotrópico)
    Id_arcilla = 1
    mapdl.mp("KXX",Id_arcilla,lambda_clay)
    mapdl.mp("KYY",Id_arcilla,lambda_clay)
    mapdl.mp("KZZ",Id_arcilla,lambda_clay)

    # Elemento finito
    Id_ladrillo = 1
    mapdl.et(Id_ladrillo,"Solid70") #element reference p.297 "SOLID70"

    # Asignar todo a ladrillo
    mapdl.vatt(Id_arcilla,0,Id_ladrillo)

    ###### MALLADO ######
    mapdl.esize(elefini)
    mapdl.vsweep(1)

    ###### CONDICIONES DE BORDE DOF CONSTRAINS ######

    Thot = 30 #C
    Tcold = 18 #C

    mapdl.asel("S", vmin=2)
    mapdl.da("all", "temp", Tcold)
    mapdl.asel("S", vmin=4)
    mapdl.da("all", "temp", Thot)
    out = mapdl.allsel()

    ###### SOLVE ###### 
    mapdl.run("/SOLU")
    mapdl.solve()
    out = mapdl.finish()

    ###### POST-PROCESSING ###### 
    mapdl.post1()
    mapdl.set("last","last")

    ###### ITERACIONES CONVECCION ###### 

    iteracion = 0

    while iteracion < 3:
        print(iteracion)
        # Calculo Tm
        a_conv = []
        for i in range(6,n_cav+2):
            a_conv = np.append(a_conv,i)
        for i in range(n_cav+3,n_cav*4+7):
            a_conv = np.append(a_conv,i)

        T_conv = np.zeros(len(a_conv))
        j = 0
        for i in  a_conv:
            mapdl.asel("S", vmin=i)
            mapdl.nsla()
            T = mapdl.post_processing.nodal_temperature() #temperaturas en cada nodo del area
            T_conv[j]= sum(T)/len(T) #temperatura promedio del area
            j = j+1

        Tm_conv = np.zeros(int(len(T_conv)/4)) #temperatura media para cada cavidad [c1,c2,c3,c4]
        j = 0
        for i in range(0,len(T_conv),4):
            Tm_conv[j] = (T_conv[i]+T_conv[i+1]+T_conv[i+2]+T_conv[i+3])/4
            j = j+1

        # Calculo parametro H (W/m2K)
        btz = 5.67e-8 #W/m2K4
        RAD = 0.83

        d = np.zeros(n_cav)
        for i in range(0,len(d)):
            d[i]= matriz[i,3] #m

        b1 = np.zeros(n_cav)
        for i in range(0,len(b1)):
            b1[i]= matriz[i,2] #m

        hr0 = 4*btz*(np.array(Tm_conv+273.15)**3)

        ha = np.zeros(n_cav)
        for i in range(0,len(ha)):
            ha[i] = max(0.025/d[i],1.25)

        hr = hr0/((1/RAD)+(1/RAD)-2+(2/(1+np.sqrt(1+(d**2)/(b1**2))-d/b1)))

        H = hr + ha #W/m2K

        # SURFACE LOAD CONVECTION

        mapdl.run("/SOLU")
        mapdl.sfcum("all","REPL") #reemplaza las cargas en areas

        j = 0
        for i in range(7,len(T_conv),4):
            mapdl.asel("S", vmin=i,vmax=i+3) #selecciona areas
            mapdl.nsla() #selecciona nodos asociados a las areas seleccionadas
            mapdl.sf("all","conv",H[j],Tm_conv[j])
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
