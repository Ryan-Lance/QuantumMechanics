from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

L = 10
T = 20

dx = 0.001
dt = 0.05

Nx = int(L/dx)
Nt = int(T/dt)

m = 1
hbar = 1
w = 1


x = np.arange(-L/2, L/2, dx)

def V(x):
	return 0.5 * w**2 * (x-L/2)**2 

def normalize(psi):
	total = 0
	for i in range(Nx):
		total += abs(psi[i])**2 * dx
	for j in range(Nx):
		psi[j] /= np.sqrt(total)
	return psi
	
def normCoeff(cn_list):
	total = 0
	N = range(len(cn_list))
	for i in N:
		total += abs(cn_list[i])**2
	for j in N:
		cn_list[j] /= np.sqrt(total)
	return cn_list

def state(n):
	psi = np.zeros(Nx)
	psi[0] = 0
	psi[1] = dx
	E = hbar*w*(n+0.45)
	dE = 1.0
	done = False
	while not done:
		psi, nodes = buildpsi(E, psi)
		if nodes <= n:
			E += dE
			#print "Enew = ", E
		if nodes > n:
			if dE > 0.0000001:
				#print "Too big"
				E -= dE
				dE /= 2.
			else:
				done = True
	return normalize(psi)

def HOeigenstate(n):
	psi = np.zeros(Nx)
	psi[0] = 0
	psi[1] = dx
	E = hbar*w*(n+0.5)
	psi, nodes = buildpsi(E, psi)
	return normalize(psi)

def buildpsi(E, psi):
        nodes = 0
        for i in np.arange(2, Nx-1):
                psi[i+1] = 2*psi[i]-2*m*(dx**2)*(E-V(i*dx))/(hbar**2)*psi[i]-psi[i-1]
                if (psi[i] > 0 and psi[i+1] < 0) or (psi[i] < 0 and psi[i+1] > 0):
                        nodes += 1
        return psi, nodes


states = []
for k in range(5):
        states.append(HOeigenstate(k))
        print "Built state ", k

# Build superposition state
coefficients = normCoeff([1,2,3,4,1])

superpositionR = np.zeros((Nt, Nx)) # psi[time][position]
superpositionI = np.zeros((Nt, Nx)) # psi[time][position]
for t in range(Nt):
        for n in range(5):
                En = hbar*w*(n + 0.5)
                first = coefficients[n]*states[n]
                superpositionR[t] += first*np.cos(-En*t*dt/hbar)
                superpositionI[t] += first*np.sin(-En*t*dt/hbar)
new = []
for i in range(5):
        new.append(round(coefficients[i],2))
plt.ion()
for t in range(Nt):
	plt.cla()
	plt.plot(x, superpositionR[t], '--', label = "Re[Psi]")
	plt.plot(x, superpositionI[t], '--', label = "Im[Psi]")
	plt.plot(x, superpositionR[t]**2+superpositionI[t]**2)
	TEXT = "Superposition: "+ str(new) + "  Time: " + str(t*dt)
	plt.text(0.1-L/2, 0.9, TEXT)
	plt.ylim([-1, 1])
	plt.xlim([-L/2, L/2])
	plt.legend()
	plt.draw()
	plt.pause(1e-3)
plt.ioff()
'''
class HOpotential():
	
	def __init__(self):
		pass
	
	def V(x):
		self.V = 0.5 * w**2 * (x-L/2)**2 
		return V
	

class wavefunction():
	
	def __init__(self, function, potential = HOpotential):
		
		self.normalize()				# normalize c_n's
		self.projectState(HOpotential)	# project state onto potential
			
	
	def normalize():
		total = 0
		for i in states:
			total += (states[i][0])**2
		for j in states:
			states[i][0] /= total
	
	def projectState():
		pass
		
'''
