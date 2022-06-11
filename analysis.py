# - - - - - - - - - - - - - - - - - - - - - - - -
# analysis.py
# 2D conductive heat transfer analysis
# - - - - - - - - - - - - - - - - - - - - - - - -
# Author: H. A. Güler
# - - - - - - - - - - - - - - - - - - - - - - - -
# Last updated: 2022-02-08
# - - - - - - - - - - - - - - - - - - - - - - - -
# Glossary
# mtl: material
# L: length
# A: area
# T: temperature
# V: volume
# m: mass
# pos: position
# heat_cpc: heat capacity (not specific)
# spec_heat: specific heat capacity
# n: number (of something)
# _i: index (of whatever comes before the underscore)
# d: near-infinitesimal amount (of something)
# - - - - - - - - - - - - - - - - - - - - - - - -
# Unless otherwise specified:
# All length units are in mm
# All area units are in mm2
# All volume units are in mm3
# All temperature units are in K
# All mass units are in g
# All time units are in s
# - - - - - - - - - - - - - - - - - - - - - - - -

import math
import matplotlib.pyplot as plt

from meshing import *
from material import SS304L, CuCrZr

SS304L = SS304L()
CuCrZr = CuCrZr()

#mesh = create_equal_cell_mesh(30, 20, 1, 30, 30, 10, CuCrZr, (273+25), (273+550), (273+400), (273+25), (273+25))
mesh = create_equal_cell_mesh_left_fixed(30, 20, 1, 30, 30, 10, CuCrZr, (273+25), (273+550))

def dist(pos1, pos2):
    return math.sqrt((pos2[0] - pos1[0])**2 + (pos2[1] - pos1[1])**2 + (pos2[2] - pos1[2])**2)

def conduct_differential_heat_on_mesh(mesh, dt):
    n_x = mesh.n_x
    n_y = mesh.n_y
    n_z = mesh.n_z

    dTs = []

    for z_i in range(len(mesh.cells)):
        for y_i in range(len(mesh.cells[0])):
            for x_i in range(len(mesh.cells[0][0])):
                cell_i_x, cell_i_y, cell_i_z = x_i, y_i, z_i
                cell = mesh.get_cell(x_i, y_i, z_i)

                if not cell.is_fixed():
                    cell_left = mesh.get_cell(cell_i_x - 1, cell_i_y, cell_i_z)
                    cell_right = mesh.get_cell(cell_i_x + 1, cell_i_y, cell_i_z)
                    cell_front = mesh.get_cell(cell_i_x, cell_i_y - 1, cell_i_z)
                    cell_rear = mesh.get_cell(cell_i_x, cell_i_y + 1, cell_i_z)

                    dx = dist(cell.get_pos(), cell_left.get_pos()) * 0.001 # convert from mm to m
                    dy = dist(cell.get_pos(), cell_front.get_pos()) * 0.001 # convert from mm to m

                    T_c = cell.get_T()
                    
                    if cell_left:
                        T_left = cell_left.get_T()
                    else:
                        T_left = T_c

                    if cell_right:
                        T_right = cell_right.get_T()
                    else:
                        T_right = T_c

                    if cell_front:
                        T_front = cell_front.get_T()
                    else:
                        T_front = T_c

                    if cell_rear:
                        T_rear = cell_rear.get_T()
                    else:
                        T_rear = T_c
                    
                    k = cell.get_thermal_conductivity()
                    rho = cell.get_density()
                    spec_heat = cell.get_spec_heat()

                    alpha = k/(rho * spec_heat) # thermal diffusivity
                    dT = alpha * (((T_left - 2*T_c + T_right)/(dx**2)) + ((T_front - 2*T_c + T_rear)/(dy**2))) * dt
                    dTs.append(dT)
                else:
                    # boundary cell, no temperature change
                    dTs.append(0)

    return dTs

def conduct_heat_in_time_interval(mesh, time_end, dt=0.001):
    time = 0
    while time < time_end:
        
        dTs = conduct_differential_heat_on_mesh(mesh, dt)

        i = 0
        for z_i in range(len(mesh.cells)):
            for y_i in range(len(mesh.cells[0])):
                for x_i in range(len(mesh.cells[0][0])):

                    cell = mesh.get_cell(x_i, y_i, z_i)
                    cell.T += dTs[i]

                    i += 1

        time += dt

conduct_heat_in_time_interval(mesh, 10, 0.01)

def plot_mesh_T(mesh):
    xs = []
    ys = []
    Ts = []

    for x_i in range(mesh.n_x):
        xs.append(mesh.get_cell(x_i, 0, 0).get_pos()[0])

    for y_i in range(mesh.n_y):
        ys.append(mesh.get_cell(0, y_i, 0).get_pos()[1])

    for y_i in range(mesh.n_y):
        Ts.append([])
        for x_i in range(mesh.n_x):
            cell = mesh.get_cell(x_i, y_i, 0)
            T = cell.get_T()
            Ts[y_i].append(T)

    clrplot = plt.contourf(xs, ys, Ts)
    cplot = plt.contour(xs, ys, Ts, colors="k")
    plt.title("Temperature Gradient (K)")
    plt.clabel(cplot, fontsize=10)
    plt.xlabel("X position (mm)")
    plt.ylabel("Y position (mm)")
    plt.show()

plot_mesh_T(mesh)
