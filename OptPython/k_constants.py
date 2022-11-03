# we distinguish between L and l or W and w in our notation
# as each one represents different things

# number of cavities per row in the brick
# each entry is the number per row
# ie: [4, 3, 4] or [3, 3, 3]
CAVITIES_PER_ROW = [4, 3, 4, 3, 4]

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

# Btz constant
BTZ = 5.67e-8 #W/m2K4

# RADiation inside the cavities
RAD = 0.83