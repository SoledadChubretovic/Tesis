import numpy as np
from k_constants import (
    CAVITIES_PER_ROW,
    L,
)

# returns the index of the first cavity of the row 
# and the index of the first cavity of the next row
# example: row = 2 has cavities [7,8,9,10]
# index_start,index_end = row_limits(2) = 7,11
def row_limits(row: int) -> tuple:
    return sum(CAVITIES_PER_ROW[:row]), sum(CAVITIES_PER_ROW[:row + 1])

# returns the index of the cavity given the row and column
# example: row = 2 has cavities [7,8,9,10]
# index_start,_ = li_index(2,1) = 8
def li_index(row: int, col: int) -> int:
    index_start, _ = row_limits(row)
    return index_start + col

# returns thickness of the partition wall in y axis
# all partitions walls have the same thickness in y direction
# because the value w of all cavities is the same 
# (w: dimension of the cavity in y axis)
# the number of partition walls in a column or row is always one more 
# than the number of cavities in that column or row
def tabique_w(W: int, w: int) -> float:
    n_cav_row = len(CAVITIES_PER_ROW)
    return (W - n_cav_row * w) / (n_cav_row + 1)
    
# returns thickness of the partition wall in x axis
# in x axis, all partitions walls in a given row have the same thikness
# value on l for each cavity is different 
# (l: dimension of cavity in x axis)
# the number of partition walls in a column or row is always one more 
# than the number of cavities in that column or row
def tabique_l(row: int, x: list) -> float:
    index_start, index_end = row_limits(row)
    return (L - sum(x[index_start: index_end])) / (CAVITIES_PER_ROW[row] + 1)

# returns the position of the center of a cavity in y axis 
def pos_y(row: int, tabique_w: float, w: int):
    return tabique_w + (w / 2) + row * (tabique_w + w)

# returns the position of the center of a cavity in x axis
def pos_x(row: int, col: int, x: list) -> float:
    index_start, index_end = row_limits(row)
    if col == 0:
        return tabique_l(row, x) + x[index_start] / 2

    pos_x_previous = pos_x(row, col - 1, x)
    index = li_index(row, col)
    cavity_previous = x[index - 1]
    cavity_current = x[index]
    return pos_x_previous + cavity_previous / 2 + cavity_current / 2 + tabique_l(row, x)

# returns (center_x, center_y, dimension_x, dimension_y) for all cavities of a brick
# ANSYS command blc5(center_x, center_y, dimension_x, dimension_y)
# blc5 command creates a rectangle given the position of the center and the dimensions in each axis.
def generate_cavitiy_ansys_parameters(x: list, W: int, w: int):
    vectors = []
    for row_idx, row_value in enumerate(CAVITIES_PER_ROW):
        tw = tabique_w(W, w)
        center_y = pos_y(row_idx, tw, w)
        for col in range(row_value):
            center_x = pos_x(row_idx, col , x)
            l_i = np.array([center_x, center_y, x[li_index(row_idx, col)], w])
            vectors.append(l_i)
    return np.array(vectors)
            
            