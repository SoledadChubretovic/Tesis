#%% matrix
import numpy as np
from functions import (
    generate_cavitiy_ansys_parameters,
    brick_corners_positions,
    plot_polygon
)
from k_constants import (
    L,
)

brick1 = [8.44662,  1.10005,  10,  43,  25,  25,  34,  34,  11,  25,  40,  25,  17,  13,  20,  23,  14,  46,  12,  298,  0.61]
brick2 = [3.46,     1.64516,  11,  42,  25,  26,  34,  31,  13,  26,  40,  14,  17, 13,  20,  23,  15,  46,  12,  140,  0.61]

x = brick2[2:]

_w_pos = len(x) - 3
_W_pos = len(x) - 2

# width of the brick
W = x[_W_pos]  # mm
# width of every cavity 
w = x[_w_pos]  # mm

b = np.array([L/2, W/2, L, W])/1000  # m
matrix = generate_cavitiy_ansys_parameters(x, W, w) / 1000  # m

draw_brick = brick_corners_positions(b, matrix)

from shapely.geometry import Polygon
import matplotlib.pyplot as plt
# Input polygon with two holes
# (remember exterior point order is ccw, holes cw else
# holes may not appear as holes.)

polygon = Polygon(shell = draw_brick[0],
                  holes = draw_brick[1:])

fig, ax = plt.subplots()
plot_polygon(ax, polygon, facecolor = 'white', edgecolor = 'red')

#%%
with open('/Users/macbookair/Desktop/Tesis-main/OptPython/1_Model/OptPython_Log.txt','r') as f:
    lines = f.read().splitlines()

