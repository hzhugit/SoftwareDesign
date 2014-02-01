# This function draws a grid.

import sys

row = int(sys.argv[1])
col = int(sys.argv[2])
length = int(sys.argv[3])

def draw_grid(row, col, length):
    for i in range(row):
        print col * ("+" + length * "-") + "+"
        for i in range(length):
            print col * ("|" + length * " ") + "|"
    print col * ("+" + length * "-") + "+"

draw_grid(row, col, length)