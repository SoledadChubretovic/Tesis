#%% matrix
from shapely.geometry import Polygon
import numpy as np
import matplotlib.pyplot as plt
import time
from functions import (
    generate_cavitiy_ansys_parameters,
    brick_corners_positions,
    plot_polygon
)
from k_constants import (
    L,
)

brick1 = [4.86769,  1.47129,  20,  66,  24,  35,  15,  35,  15,  24,  22,  10,  30,  33,  15,  26,  66,  30,  18,  13,  198,  0.59]
brick2 = [7.16879,  1.32368,  20,  67,  24,  35,  15,  39,  15,  25,  35,  13,  29,  34,  15,  26,  66,  30,  18,  16,  284,  0.59]

x = brick1[2:]

_w_pos = len(x) - 3
_W_pos = len(x) - 2

# width of the brick
W = x[_W_pos]  # mm
# width of every cavity 
w = x[_w_pos]  # mm

b = np.array([L/2, W/2, L, W])/1000  # m
matrix = generate_cavitiy_ansys_parameters(x, W, w) / 1000  # m

draw_brick = brick_corners_positions(b, matrix)

# Input polygon with two holes
# (remember exterior point order is ccw, holes cw else
# holes may not appear as holes.)

polygon = Polygon(shell = draw_brick[0],
                  holes = draw_brick[1:])

fig, ax = plt.subplots()
plot_polygon(ax, polygon, facecolor = 'white', edgecolor = 'red')
#%%
#plt.savefig(r"C:\Users\nchubretovic\OneDrive - Entel\Escritorio\Sole\Tesis\1_Model\bricks_images\Brick_" + time.strftime("%Y-%m-%d %H.%M.%S") + ".png")



# %%
