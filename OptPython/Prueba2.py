#%% ACTIVAR LICENCIA Y LIBRERIAS

# Unidades: m,W,C

import numpy as np
#from ansys.mapdl.core import launch_mapdl
#mapdl = launch_mapdl()

# %% CALCULO TERMICO LADRILLO
import numpy as np
mapdl.finish()
mapdl.clear()
mapdl.prep7()

###### INPUTS DESDE QUADRATICFUNCTION.M ######

# [centro x, centro y, largo x, largo y]

#### GEOMETRY AND MATERIAL VALUES ####

x =[16       , 154      ,   66   ,     104      ,   32     ,    26    ,     37 ,        44]

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
dimz = 0.02 #m
# H = 100 #mm
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
p_brick = p_clay*An/(W*L) #MPa
p_brickG2 = 15 #Mpa

# Thermal parameters
k = 0.15 #W/mK
U_Santiago = 5 #W/m2K

# [centro x, centro y, largo x, largo y]
b1 = np.array([L/2,W/2,L,W])/1000 #m
c1 = np.array([tl12+l1/2,tw13+w1/2,l1,w1])/1000 #m
c2 = np.array([2*tl12+l1+l2/2,tw24+w2/2,l2,w2])/1000 #m
c3 = np.array([tl34+l3/2,2*tw13+w1+w3/2,l3,w3])/1000 #m
c4 = np.array([2*tl34+l3+l4/2,2*tw24+w2+w4/2,l4,w4])/1000; #m


matriz = np.array([c1,c2,c3,c4])
#dimz = 0.02 #m
#mapdl.aplot()
k = 0.5 #W/mC

###### GEOMETRIA ######

brick_perimeter = mapdl.blc5(b1[0],b1[1],b1[2],b1[3])
for row in matriz:
    mapdl.blc5(row[0],row[1],row[2],row[3])

cavities = mapdl.asba(brick_perimeter, 'all')
mapdl.allsel()
mapdl.vext(cavities,dz = dimz) #m
#mapdl.vplot()
mapdl.aplot(show_area_numbering=True)

###### ASIGNACION DE ATRIBUTOS ###### 

# Material (Conductividad termica, material isotrópico)
Id_arcilla = 1
mapdl.mp("KXX",Id_arcilla,k)
mapdl.mp("KYY",Id_arcilla,k)
mapdl.mp("KZZ",Id_arcilla,k)

# Elemento finito
Id_ladrillo = 1
mapdl.et(Id_ladrillo,"Solid70") #element reference p.297 "SOLID70"

# Asignar todo a ladrillo
mapdl.vatt(Id_arcilla,0,Id_ladrillo)

###### MALLADO ######

mapdl.vsweep(1)
#mapdl.eplot()

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
mapdl.post_processing.plot_nodal_temperature()

###### ITERACIONES CONVECCION ###### 

iteracion = 0

while iteracion < 7:
    # Calculo Tm
    # Necesito temperaturas en area 8,10,12,14,16,18,20,22
    a_conv = [8,10,12,14,16,18,20,22] #que coincida con el orden de la geometria
    T_conv = np.zeros(len(a_conv))
    j = 0

    for i in a_conv:
        mapdl.asel("S", vmin=i)
        mapdl.nsla()
        T = mapdl.post_processing.nodal_temperature() #temperaturas en cada nodo del area
        T_conv[j]= sum(T)/len(T) #temperatura promedio del area
        j = j+1

    Tm_conv = np.zeros(int(len(T_conv)/2)) #temperatura media para cada cavidad [c1,c2,c3,c4]
    j = 0

    for i in range(0,len(T_conv),2):
        Tm_conv[j] = (T_conv[i]+T_conv[i+1])/2
        j = j+1

    # Calculo parametro H (W/m2K)
    btz = 5.67e-8 #W/m2K4
    RAD = 0.83

    d = np.zeros(4)
    for i in range(0,len(d)):
        d[i]= matriz[i,3] #m

    b = np.zeros(4)
    for i in range(0,len(b)):
        b[i]= matriz[i,2] #m

    hr0 = 4*btz*(np.array(Tm_conv+273.15)**3)

    ha = np.zeros(4)
    for i in range(0,len(ha)):
        ha[i] = max(0.025/d[i],1.25)

    hr = hr0/((1/RAD)+(1/RAD)-2+(2/(1+np.sqrt(1+(d**2)/(b**2))-d/b)))

    H = hr + ha #W/m2K
    # print(iteracion)
    # print(H)
    # print(Tm_conv)

    # SURFACE LOAD CONVECTION

    mapdl.run("/SOLU")
    mapdl.sfcum("all","REPL") #reemplaza las cargas en areas

    #cavidad 1
    mapdl.asel("S", vmin=7,vmax=10) #selecciona areas
    mapdl.nsla() #selecciona nodos asociados a las areas seleccionadas
    mapdl.sf("all","conv",H[0],Tm_conv[0])

    #cavidad 2
    mapdl.asel("S", vmin=11,vmax=14) #selecciona areas
    mapdl.nsla() #selecciona nodos asociados a las areas seleccionadas
    mapdl.sf("all","conv",H[1],Tm_conv[1])

    #cavidad 3
    mapdl.asel("S", vmin=15,vmax=18) #selecciona areas
    mapdl.nsla() #selecciona nodos asociados a las areas seleccionadas
    mapdl.sf("all","conv",H[2],Tm_conv[2])

    #cavidad 4
    mapdl.asel("S", vmin=19,vmax=22) #selecciona areas
    mapdl.nsla() #selecciona nodos asociados a las areas seleccionadas
    mapdl.sf("all","conv",H[3],Tm_conv[3])

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
flux_tot = flux_area*(dimz*b1[2])#W

###### CONDUCTIVIDAD EQUIVALENTE ###### 

keq = flux_area*b1[3]/(Thot-Tcold) #W/mK
print(keq)
# ###### TRANSMITANCIA DEL MURO ######

RsiRse = 0.17 #m2K/W
k_mort_pega = 0.23 #W/mK
prop_ladrillo = 0.84
prop_mort = 1-prop_ladrillo
C_muro = (prop_ladrillo*keq+prop_mort*k_mort_pega)/b1[3] #W/m2K
U_muro = 1/(RsiRse+1/C_muro) #W/m2K
print(U_muro)

# %% LIBERAR LICENCIA - STOP MAPDL

mapdl.exit()

# %%
