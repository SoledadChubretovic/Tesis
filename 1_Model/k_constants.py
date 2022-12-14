#%% CONSTANTS
# we distinguish between L and l or W and w in our notation
# as each one represents different things

# number of cavities per row in the brick
# each entry is the number per row
# ie: [4, 3, 4] or [3, 3, 3, 3, 3]
CAVITIES_PER_ROW = [3,4,3,4,3]

# total amount of cavities in the brick
N_CAVITIES = sum(CAVITIES_PER_ROW)

# total amount of variables
# this is defined as sum(li) + w + W + lambda_clay
N_VARIABLES = N_CAVITIES + 3

# length of the brick (x axis) in mm
L = 154 # mm

# thickness limits for the brick
# the constants are defined in terms of half of the thickness
# to make easier future calculations in terms of how it is modeled
# W2_MIN changes if we add more rows
W2_MIN = 200 / 10 # mm  70 * 2 = 140
W2_MAX = 300 / 10 # mm 150 * 2 = 300

# high of the brick (z axis) in mm
H = 113 # mm
DIM_Z = H / 1000 # m

# minimum length of partition walls in mm
T_MIN = 10 # mm

# clay thermal transmitance limits in W/mK
# the number is *100 to use discrete calculations and simplify
# the optimization process
LAMBDA_CLAY100_MIN = 60 #50 # W/mK
LAMBDA_CLAY100_MAX = 60 #62 # W/mK


# maximum width for cavities in mm
# setting all the other cavities and partition walls to t_min (minimum value)
w_MAX = ((W2_MAX * 2) - (len(CAVITIES_PER_ROW) + 1) * T_MIN) / len(CAVITIES_PER_ROW) # mm

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
# 0.0055 is the better choice for time and precision
ELEMENT_SIZE = 0.007

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
# based on result observation
# not the best way to end the optimization but I could not make tolerance termination work
N_GENERATION = N_VARIABLES * 5

# population size
POPULATION_SIZE = 60

# number of offsprings
# for now, set equal to POPULATION_SIZE
N_OFFSPRINGS = POPULATION_SIZE

# big number to assign to invalid geometries
U_MURO_INVALID = 50

# %%
