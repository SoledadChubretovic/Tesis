#%% otro
import numpy as np
from Termico4 import ThermalAnalysis

x = np.ones(17)*10

#### GEOMETRY AND MATERIAL VALUES ####

L = 320 #mm (Hacer calzar con L)
H = 113 #mm
W = 154
t_min = 10 #mm ancho minimo de tabiques (Hacer calzar con xl)
n_cav_L = 4
n_cav_W = 4
n_cav = n_cav_L*n_cav_W

# Areas
Ah = x[-3]*(sum(x)-W-x[-3])
Ab = W*L
An = abs(Ab-Ah) #mm2 porque en casos de que Ahuecos mayor que Abruta An era negativo

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

# Thermal parameters
lambda_clay = 0.5 #W/mK
U_Santiago = 5 #W/m2K

# Structural parameters
rho_clay = 2.0841e-6*lambda_clay + 0.55042e-6 #1.6516e-6 kg/mm3
resist_clay = 76.766*lambda_clay-26.362 #MPa
resist_brick = resist_clay*An/Ab #MPa
resist_brickG2 = 10 #Mpa

# [centro x, centro y, largo x, largo y]
b = np.array([L/2,W/2,L,W])/1000 #m
tmin = t_min/1000 #m
dimz = H/1000 #m

b = np.array([L/2,W/2,L,W])/1000 #m
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

suma = sum(x[0:4])
print(x[0],x[1],x[2],x[3],suma)

# %%
