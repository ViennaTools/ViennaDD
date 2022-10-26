
from PSSolver import PSSolver2D
from PSSolver import PSSolver3D
import matplotlib.pyplot as plt
import math
from timeit import default_timer as timer
import numpy as np

"""
Solves the Poisson equation for a known solution in 2D. The maximum occurring error is then printed.
"""


if __name__ == '__main__':
    # u = sin(x)*cos(y)
    # uxx + uyy = -2sin(x)*sin(y)
    nx = 20
    ny = 20
    meshgrid_x = np.linspace(0,np.pi, nx)
    meshgrid_y = np.linspace(0,np.pi, ny)

    solver = PSSolver2D()
    xa, xb = np.meshgrid(meshgrid_x, meshgrid_y)

    meshgrid = []
    for i,j in zip(xa, xb):
       for xkoord, ykoord in zip(i,j): 
           meshgrid.append([xkoord, ykoord])
    solver.setSimulationGrid(meshgrid)
    cells = []
    for i in range(ny-1):
        for j in range(nx-1):
            cells.append([ny*i+j, ny*i+j+1, ny*(i+1)+j+1, ny*(i+1)+j])
    
    simulationGrid = []
    for i in range(len(cells)):
        simulationGrid.append(i)

    solver.addSemiconductorSegment(simulationGrid, "simulationDomain")
    solver.import_rectangular_mesh(cells, nx, ny)

    f = []
    for i in meshgrid:
        f.append(2*np.sin(i[0])*np.cos(i[1]))


    kathode_points = []
    for i in range(ny):
        kathode_points.append(i*nx)
    
    anode_points = []
    for i in range(ny):
        anode_points.append(i*nx+nx-1)  
    
    solver.setChargeDensity(f)
    solver.addContactPoints(anode_points, "Anode")
    solver.addContactPoints(kathode_points, "Katode")

    solver.applyVoltageSegment(0, "Katode")
    solver.applyVoltageSegment(0, "Anode")
    
    solver.test()
    
    result = solver.solvePoisson()

    correct =  []
    for i in meshgrid:
        correct.append(np.sin(i[0])*np.cos(i[1]))

    max_error = 0
    for i in range(len(correct)):
        if abs(correct[i]-result[i]) > max_error:
            max_error = abs(result[i]-correct[i])
    print("max error ", max_error)    