#Wyatt Davis
#Final Project (M.N. Computational Physics, problem 9.8)

#Import all relevant functions.
import numpy as np  
import math as math
import cmath as cmath
from banded import banded
import matplotlib.pyplot as plt




#Prompt user to input values pretaining to the data collection
num=input("Input the desired number of trial runs: ")
E=input("Input a seed value for the energy (Ev) for the free particle: ")
Int=input("Input the desired interval of trial energies: ")
For=input("Input the time for each trial (atto seconds): ")


#Assign a variable to the length of the well. 
L=1e-8
#Assign a variable to the sigma value.
sig=1e-10
#Assign a variable to the k value. 
k=5e10
#Assign a variable to the initial x position.
xi=0.
#Assign a variable to the length of the well.
xf=L
#Assign a variable to number of intervals.
N=3000
#Assign a variable to the grid spacing.
a=(xf-xi)/N
#Assign a variable to the time step.
h=1e-18
#Assign a value to the initial x position.
x0=L/2
#Assign a value to hbar. 
h_b=1.053e-34
#Assign a value to the mass of the electron.
m=9.109e-31

# This function calculates the t=0 array for psi for all x between 0 and L. Keeps x=0,L at zero in order to satisfy the boundary conditions. 
def psi_i(xi,xf,x0,sig,k,N):
	a=(xf-xi)/N
	psi=np.zeros([N+1,1], complex)
	for i in range(0,N):
		psi[i,0]=math.exp(-(((i*a)-x0)**2)/(2*(sig**2)))*cmath.exp(1j*k*i*a)
	return psi

#Define the position function that I will utilize in updating the position of the curve object. 
def pos(phi,x):
    N=len(phi)
    func=np.zeros([N,3],float)
    func[:,0]=x[:,0]
    func[:,1]=abs(phi[:,0])*3*(1e-9)
    return func

#Define the function that will calculate the relevant 'v' vector argument for the banded method. 
def v(phi,b1,b2):
    v=np.zeros([N+1,1],complex)
    v[0,0]=phi[0,0]*b1+b2*phi[1,0]
    for i in range(1,N-1):
        v[i,0]=b1*phi[i,0]+b2*(phi[i+1,0]+phi[i-1,0])
    v[-1,0]=b1*phi[-1,0]
    return v 

#Function that returns the most probable x position of the particle
def ClassicalPos(phi,x):
    EL=np.argmax(abs(phi[:,0]))
    pos=x[EL,0]
    return pos
    

#Initialize the array.
A=np.zeros([3,N+1],complex)

#Define the values to go into the banded matricies

a1=1+h*(1j*h_b/2/m/(a**2))
a2=-h*(1j*h_b/4/m/(a**2))
b1=1-h*(1j*h_b/2/m/(a**2))
b2=h*(1j*h_b/4/m/(a**2))

#Filling the A-array with constants, structure condusive to banded function. 
A[0,:],A[1,:],A[2,:]=a2,a1,a2

#Initialize x-position array
x=np.zeros([N+1,1],float)

#Fill x-position array. 
for il in range(0,N+1):
    x[il,0]=(-L/2)+il*a

#Getting inputted energies in terms of joules
E*=1.602176e-19
Int*=1.602176e-19

Kx=[]
V=[]
'''
Loop that takes multiple linear regressions of classical position v. time data and
appends best fit momentums and cooresponding k values to their respective lists.
'''
for n in range(num):
        pos=[]
        t=[]
        time=0
        k=math.sqrt(2*m*E/h_b**2)
        phi=psi_i(xi,xf,x0,sig,k,N)
        for i in range(For):
                phi=banded(A,v(phi,b1,b2),1,1)
                pos.append(ClassicalPos(phi,x))
                t.append(time)
                time+=h
        vel,b=np.polyfit(t,pos,1)
        Kx.append(k)
        V.append(vel*m)
        E+=Int
'''
Finds linear regression of momentum v. k data. According to the De-Broglie momentum relation: m*v=k*hbar,
the slope of this graph should be hbar.
'''
Fit=np.polyfit(Kx,V,1)

#Output the percent of hbar obtained.
Percent=(Fit[0]/h_b)*100
Per=str(Percent)
print("The experimentally obtained value for hbar is", Per, "percent of the accepted value.")  
plt.plot(Kx,V)
plt.show()

### NEXT STEP ###

#How to calculate error?



