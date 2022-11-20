#%% CONSTANTS
# we distinguish between L and l or W and w in our notation
# as each one represents different things

# number of cavities per row in the brick
# each entry is the number per row
# ie: [4, 3, 4] or [3, 3, 3, 3, 3]
CAVITIES_PER_ROW = [4, 4, 4, 4]

# total amount of cavities in the brick
N_CAVITIES = sum(CAVITIES_PER_ROW)

# total amount of variables
# this is defined as sum(li) + w + W + lambda_clay
N_VARIABLES = N_CAVITIES + 3 

# length of the brick (x axis) in mm
L = 154 # mm

# high of the brick (z axis) in mm
H = 113 # mm
DIM_Z = H / 1000 # m

# minimum length of partition walls in mm
T_MIN = 10 # mm

# thickness limits for the brick
# the constants are defined in terms of half of the thickness
# to make easier future calculations in terms of how it is modeled
W2_MIN = 140 / 2 # mm  70 * 2 = 140
W2_MAX = 300 / 2 # mm 150 * 2 = 300

# clay thermal transmitance limits in W/mK
# the number is *100 to use discrete calculations and simplify
# the optimization process
LAMBDA_CLAY100_MIN = 50 # W/mK
LAMBDA_CLAY100_MAX = 62 # W/mK

# maximum lenght for cavities in mm
# setting all the other cavities and internal walls to t_min (minimum value)
l_MAX = L - 2 * min(CAVITIES_PER_ROW) * T_MIN # mm

# maximum width for cavities in mm
# setting all the other cavities and internal walls to t_min (minimum value)
w_MAX = W2_MAX - 2 * len(CAVITIES_PER_ROW) * T_MIN # mm

# Stefan-Boltzmann constant
BTZ = 5.67e-8 #W/m2K4

# Radiation fator inside the cavities: emissivity of surfaces
# air radiation emissivity form literature
RAD = 0.83

# Radiation factor for inner and outer faces of a wall
# based NCh853
RsiRse = 0.17  # m2K/W

# thermal conductivity of glue mortar (experimetal value) in W/mK
λ_GLUE_MORTAR = 0.23  # W/mK

# brick proportion in 1 square meter wall
# this value is considered with 12 mm glue mortar between bricks
# this value is sonsiderad for the given L and H
PROP_BRICK = 0.84

# glue mortar proportion in 1 quare meter wall
prop_mortar = 1 - PROP_BRICK

# cold area index
# area of the brick that 18 celcius degrees temperature is applied
COLD_AREA_INDEX = 2

# hot area index
# area of the brick that 30 celcius degrees temperature is applied
HOT_AREA_INDEX = 4

# temperatures applied to generate heat flux
T_HOT = 30 # degrees celcius
T_COLD = 18 # degrees celcius

# finite element size
ELEMENT_SIZE = 0.0055

# number of constraints
N_CONSTRAINTS = 3

# number of objectives
N_OBJECTIVES = 2

# normative compression resistance
# in chile bricks have to be grade G2 (11 MPa minimum)
# most bricks are grade G1 (15 MPa minimum) so we set use this value
# See NCh1928 and NCh169
G2 = 15

# number of generation for termination criterion
N_GENERATION = 30

# population size
POPULATION_SIZE = 50

# number of offsprings
N_OFFSPRINGS = 50

# big number to assign to invalid geometries
U_MURO_INVALID = 50

# %%
