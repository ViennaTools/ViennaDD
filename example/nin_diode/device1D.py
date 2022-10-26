import numpy as np
import matplotlib.pyplot as plt
from PSSolver import PSSolver1D 
import time


"""
simulates a nin-diode in 1D. No recombination is used. As a semiconductor material silicon is used
"""


donor = []
aceptor = []
plotx = []
plotpot= []
plotelec = []
plothole = []
plotcharge = []
plotdonor = [] 
plotaceptor = []
plotfield = []

plotquasielec = []
plotquasihole = []



T = 300
q = 1.602e-19
e = 8.854e-12*11.9

Na = 1e18
Nd = 1e18
kb = 1.381e-23
Ec = 1.12/2
Ev = -1.12/2
Nc = 2.8e19
Nv = 1.04e19
Tau_n0  = 5e-7
Tau_p0 = 5e-7  
Nsrh_n = 5e16
Nsrh_p =  5e16
mu_n = 500e-4
mu_p = 50e-4
apply = 0.8


x_max = 6e-7
number_of_points = 400
mesh = np.linspace(0,x_max, number_of_points)

for x in mesh :
    if x <  2e-7 or x > 4e-7:
        aceptor.append(1e1)
        donor.append(1e19)
    else:
        aceptor.append(1e3)
        donor.append(1e17)

cells = []
for i in range(0,len(mesh)-1):
    cells.append([i, i+1])
semicCells = [*range(0,len(mesh))]
mesh = np.reshape(mesh,(-1,1))

device = PSSolver1D()
device.addSemiconductorSegment(semicCells, "Semic")
device.addContactPoints([0],"Kathode")
device.addContactPoints([len(mesh)-1],"Anode")
device.setSimulationGrid(mesh)
device.setCellGrid(cells)
device.setDonorConcentration(donor)
device.setAcceptorConcentration(aceptor)


device.setPermitivity(11.7)
device.setNumberNewtonIterations(10)
device.setNumberGummelIterations(50)
#steps*Ut is the stepwide for the stepwise applied voltage
device.setSteps(0.5)
device.setConductionBandEnergy(Ec)
device.setValenceBandEnergy(Ev)
device.setHolesDensityofStates(Nv)
device.setElecDensityofStates(Nc)
device.setTrapEnergy(1.12/2)
device.setElectronMobility(400)
device.setHoleMobility(50)
device.setSpond(False)
device.setAuger(False)
device.setSHR(False)

device.applyVoltageSegment(1, "Kathode")
device.applyVoltageSegment(0, "Anode")


device.solveDDQuasiFermi()


charge = device.getChargeDensity()
elec = device.getElectronConcentration()
hole = device.getHoleConcentration()
potential = device.getPotential()

device.test()
num = 0
for x in mesh:
    plotx.append(mesh[num])
    plotpot.append(potential[num])
    plotcharge.append(charge[num])
    plotelec.append(elec[num])
    plothole.append(hole[num])
    plotaceptor.append(aceptor[num])
    plotdonor.append(donor[num])
    num = num+1
 

plt.subplot(4,1,1)
plt.plot(plotx, plotdonor, color='b', label='donor concentration')
plt.plot(plotx, plotaceptor, color='r', label='aceptor concentration')
plt.xlabel('x in μm')
plt.ylabel('dopingin cm^-3')
plt.legend()

plt.subplot(4,1,2)
plt.plot(mesh, plotelec, color='b', label='elec')
plt.plot(mesh, plothole, color='r', label='hole')
plt.xlabel('x in μm')
plt.ylabel('densities in cm^-3')
plt.legend()

plt.subplot(4,1,3)
plt.plot(plotx, plotpot, color='b', label='built in potential')
plt.xlabel('x in μm')
plt.ylabel('Potential in V')
plt.legend()

plt.subplot(4,1,4)
plt.plot(plotx, plotcharge, color='b', label='charge density')
plt.xlabel('x in μm')
plt.ylabel('Charge density in As/m^3')
plt.legend()

plt.legend()
plt.show()

