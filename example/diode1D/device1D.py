import numpy as np
import matplotlib.pyplot as plt
from PSSolver import PSSolver1D 
import time

"""
simulates a pn-diode in 1D. Shockley-Hall-Read, spontanous and Auger recombination are used. As a semiconductor material silicon is used
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
Ec = 1.12
Ev = 0
Nc = 2.8e19
Nv = 1.04e19
Tau_n0  = 5e-7
Tau_p0 = 5e-7  
Nsrh_n = 5e16
Nsrh_p =  5e16
mu_n = 500e-4
mu_p = 50e-4
apply = 0.8


Tau_n = Tau_n0/(1+(Na+Nd)/(Nsrh_n))   
Tau_p = Tau_n0/(1+(Na+Nd)/(Nsrh_p))   

Ut = (kb*T)/q
ni = np.sqrt(Nc*Nv*np.exp((Ev-Ec)/Ut))
print("Thermal voltage ", Ut)
print("Intrinsic concentration ", ni)
Ud = Ut*np.log((Na*Nd)/(ni*ni))
print("Ud ", Ud)
N0 = (Na*Nd)/(Na+Nd)
W = np.sqrt((2*e*Ud)/(q*N0*1e6))
W_n = (W*Na)/(Nd+Na)
W_p = (W*Nd)/(Nd+Na)
print("Wn " , W_n)
print("Wp " , W_p)

 
Ldiff = np.sqrt(mu_n*Tau_n)
if Ldiff < np.sqrt(mu_p*Tau_p) :
    Ldiff = np.sqrt(mu_p*Tau_p)

print("diffusion length ", Ldiff)

x_max = 0

if x_max < W_n :
    x_max = W_n

if x_max < W_p :
    x_max = W_p

if x_max < Ldiff:
    x_max = Ldiff/10

print("depletion region ", W)
print("simulation area ", x_max)

Ldn = np.sqrt((e*Ut)/(q*Nd*1e6))
Ldp = np.sqrt((e*Ut)/(q*Na*1e6))



dx = Ldn
if dx > Ldp : 
    dx = Ldp
print("Mesh size: ", dx)
number = int(x_max/dx) + 1
print("number of points ", number)


mesh = np.linspace(0,x_max, number)

for x in mesh :
    if x <  x_max/2:
        aceptor.append(Na)
        donor.append((ni*ni)/Na)
    else:
        aceptor.append((ni*ni)/Nd)
        donor.append(Nd)

cells = []
for i in range(0,len(mesh)-1):
    cells.append([i, i+1])
semicCells = [*range(2,len(mesh)-3)]
mesh = np.reshape(mesh,(-1,1))



device = PSSolver1D()
device.addSemiconductorSegment(semicCells, "Semic")
device.addContactSegment([0,1], "Kathode")
device.addContactSegment([len(mesh)-2, len(mesh)-3], "Anode")
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
device.setElecLifetime(Tau_n)
device.setHoleLifetime(Tau_p)
device.setTrapEnergy(1.12/2)
device.setElectronMobility(400)
device.setHoleMobility(50)
device.setSpond(True)
device.setAuger(True)
device.setSHR(True)

device.applyVoltageSegment(0.3, "Kathode")
device.applyVoltageSegment(0, "Anode")

device.solveDDClassic()
charge = device.getChargeDensity()
elec = device.getElectronConcentration()
hole = device.getHoleConcentration()
potential = device.getPotential()


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

