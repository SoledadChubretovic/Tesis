#%% ACTIVAR LICENCIA Y LIBRERIAS

# Unidades: m,W,C
#from ansys.mapdl.core import launch_mapdl
#mapdl = launch_mapdl()

# %% CALCULO TERMICO LADRILLO
import numpy as np
import time

start = time.time()


mapdl.finish()
mapdl.clear()
mapdl.prep7()

###### INPUTS DESDE QUADRATICFUNCTION.M ######

# [centro x, centro y, largo x, largo y]

#### GEOMETRY AND MATERIAL VALUES ####
# Number of variables, constrains and objective functions
n_cav_L = 4
n_cav_W = 4
n_cav = n_cav_L*n_cav_W
n_var = n_cav + 3 #number of variables (li,w,W,lambda_clay)

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



x = np.ones(17)*20
x = np.append(x,80)
x = np.append(x,1)

#### GEOMETRY AND MATERIAL VALUES ####
L = 320 #mm (Hacer calzar con L)
H = 10 #mm
t_min = 10 #mm ancho minimo de tabiques (Hacer calzar con xl)
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

c1a = np.array([tla+x[0]/2 , W-(tw+x[-3]/2) , x[0] , x[-3]])/1000 #m
c2a = np.array([2*tla+x[0]+x[1]/2 , W-(tw+x[-3]/2) , x[1] , x[-3]])/1000 #m
c3a = np.array([3*tla+x[0]+x[1]+x[2]/2 , W-(tw+x[-3]/2) , x[2] , x[-3]])/1000 #m
c4a = np.array([4*tla+x[0]+x[1]+x[2]+x[3]/2 , W-(tw+x[-3]/2) , x[3] , x[-3]])/1000 #m

c1b = np.array([tlb+x[0]/2 , W-(2*tw+x[-3]*3/2) , x[4] , x[-3]])/1000 #m
c2b = np.array([2*tlb+x[0]+x[1]/2 , W-(2*tw+x[-3]*3/2) , x[5] , x[-3]])/1000 #m
c3b = np.array([3*tlb+x[0]+x[1]+x[2]/2 , W-(2*tw+x[-3]*3/2) , x[6] , x[-3]])/1000 #m
c4b = np.array([4*tlb+x[0]+x[1]+x[2]+x[3]/2 , W-(2*tw+x[-3]*3/2) , x[7] , x[-3]])/1000 #m

c1c = np.array([tlc+x[0]/2 , 2*tw+x[-3]*3/2 , x[8] , x[-3]])/1000 #m
c2c = np.array([2*tlc+x[0]+x[1]/2 , 2*tw+x[-3]*3/2 , x[9] , x[-3]])/1000 #m
c3c = np.array([3*tlc+x[0]+x[1]+x[2]/2 , 2*tw+x[-3]*3/2 , x[10] , x[-3]])/1000 #m
c4c = np.array([4*tlc+x[0]+x[1]+x[2]+x[3]/2 ,2*tw+x[-3]*3/2 , x[11] , x[-3]])/1000 #m

c1d = np.array([tld+x[0]/2 , tw+x[-3]/2 , x[12] , x[-3]])/1000 #m
c2d = np.array([2*tld+x[0]+x[1]/2 , tw+x[-3]/2 , x[13] , x[-3]])/1000 #m
c3d = np.array([3*tld+x[0]+x[1]+x[2]/2 , tw+x[-3]/2 , x[14] , x[-3]])/1000 #m
c4d = np.array([4*tld+x[0]+x[1]+x[2]+x[3]/2 , tw+x[-3]/2 , x[15] , x[-3]])/1000 #m

ca = np.concatenate((c1a,c2a,c3a,c4a))
cb = np.concatenate((c1b,c2b,c3b,c4b))
cc = np.concatenate((c1c,c2c,c3c,c4c))
cd = np.concatenate((c1d,c2d,c3d,c4d))

# Thermal parameters
lambda_clay = x[-1]/100 #W/mK

matriz = np.concatenate([c1a,c2a,c3a,c4a,c1b,c2b,c3b,c4b,c1c,c2c,c3c,c4c,c1d,c2d,c3d,c4d])

###### GEOMETRIA ######
brick_perimeter = mapdl.blc5(b[0],b[1],b[2],b[3])
for row in matriz:
    mapdl.blc5(row[0],row[1],row[2],row[3])

cavities = mapdl.asba(brick_perimeter, 'all')
mapdl.allsel()
mapdl.aplot(show_area_numbering=True)
mapdl.vext(cavities,dz = dimz) #m
mapdl.aplot(show_area_numbering=True)



end = time.time()

print("Time: ", (end-start), "s")

# %%
