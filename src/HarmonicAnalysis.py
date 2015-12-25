
#Wyatt Davis
#Final Project (M.N. Computational Physics, problem 9.8)

#Import all relevant functions.
import numpy as np  
import math as math
import cmath as cmath
from banded import banded
import matplotlib.pyplot as plt




#Prompt user to input values pretaining to the data collection

E=input("Input a seed value for the energy (Ev) for the free particle: ")
Int=input("Input the desired interval of trial energies: ")

E*=1.602176e-19

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
def v(phi,bv1,bv2):
    v=np.zeros([N+1,1],complex)
    v[0,0]=phi[0,0]*bv1[0,0]+bv2[1,0]*phi[1,0]
    for i in range(1,N):
        v[i,0]=bv1[i,0]*phi[i,0]+(bv2[i+1,0]*phi[i+1,0])+(bv2[i-1,0]*phi[i-1,0])                                     
    v[-1,0]=bv1[-1,0]*phi[-1,0]+bv2[N-1,0]*phi[N-1,0]
    return v 


#Function that returns the most probable x position of the particle
def ClassicalPos(phi,x):
    EL=np.argmax(abs(phi[:,0]))
    pos=x[EL,0]
    return pos
    

#Initialize the array.
A=np.zeros([3,N+1],complex)

#Define the values to go into the banded matricies



#Filling the A-array with constants, structure condusive to banded function. 


#Initialize x-position array
x=np.zeros([N+1,1],float)

#Fill x-position array. 
for il in range(0,N+1):
    x[il,0]=(-L/2)+il*a

def harmonicoscillator(x,k):
        L=len(x[:,0])
        V=[]
        for i in range(L):
                V.append(.5*k*(x[i,0]**2))
        return V

av1=np.zeros([N+1,1], complex)
av2=np.zeros([N+1,1], complex)
bv1=np.zeros([N+1,1], complex)
bv2=np.zeros([N+1,1], complex)
K=70
V=harmonicoscillator(x,K)
for i in range(0,N+1):
    av1[i,0]=(h/2)*((2/h)+(1j*h_b/(m*(a**2)))+(1j*V[i]/h_b))
    av2[i,0]=-1j*h*h_b/4/m/(a**2)
    bv1[i,0]=-(h/2)*((1j*h_b/m/(a**2))+(1j*V[i]/h_b)-(2/h))
    bv2[i,0]=1j*h*h_b/4/m/(a**2)

A[0,:],A[1,:],A[2,:]=av2[:,0],av1[:,0],av2[:,0]

phi=psi_i(xi,xf,x0,sig,k,N)


'''
#Getting inputted energies in terms of joules
E*=1.602176e-19
Int*=1.602176e-19

Kx=[]
V=[]
'''
time=0
t=[]
Pos=[]
while (time<2000*h):
    phi=banded(A,v(phi,bv1,bv2),1,1)
    Pos.append(ClassicalPos(phi,x))
    t.append(time)
    time+=h
plt.plot(t,Pos)
plt.show()
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
'''
Finds linear regression of momentum v. k data. According to the De-Broglie momentum relation: m*v=k*hbar,
the slope of this graph should be hbar.
'''
#Fit=np.polyfit(Kx,V,1)

#Output the percent of hbar obtained.

### NEXT STEP ###

#How to calculate error?


