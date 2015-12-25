#Wyatt Davis
#Senior Project: Harmonic Oscillator

#Import all relevant functions.
import numpy as np  
import math as math
import cmath as cmath
from visual import curve,rate,color,arrow,text,display
from banded import banded

#set scene 
scene=display(title='Wave Packet Evolution: Quantum Harmonic Oscillator', background=color.white)

#prompt user to input energy value for the particle and the spring constant for the harmonic oscillator.
E=input("Input a value for the initial kinetic energy (ev) for the particle: ")
K=input("Input a value for the spring constant of the harmonic oscillator: ")


E*=1.602176e-19
#Assign a variable to the length of the well. 
L=2.5e-9
#Assign a variable to the sigma value.
sig=1e-10
#Assign a variable to the initial x position.
xi=0
#Assign a variable to the length of the well.
xf=L
#Assign a variable to number of intervals.
N=250
#Assign a variable to the grid spacing.
a=(xf-xi)/N
#Assign a variable to the time step.
h=1e-18
#Assign a value to the initial x position.
x0=L*(3./5.)
#Assign a value to hbar. 
h_b=1.053e-34
#Assign a value to the mass of the electron.
m=9.109e-31
'''
Assign a variable to the k value for wave packet representation of particle,
utilizing classical relation between classical momentum and energy in addition
to the de Broglie Momentum relation.
'''
k=math.sqrt(2*m*E/(h_b**2))

# This function calculates the t=0 array for psi for all x between 0 and L. Keeps x=0,L at zero in order to satisfy the boundary conditions. 
def psi_i(xi,xf,x0,sig,k,N):
	a=(xf-xi)/N
	psi=np.zeros([N+1,1], complex)
	for i in range(0,N+1):
		psi[i,0]=(math.exp(-(((i*a)-x0)**2)/(2*(sig**2)))*cmath.exp(1j*k*i*a))
	return psi

#Define the position function that I will utilize in updating the position of the curve object. 
def pos(phi,x):
    N=len(phi)
    func=np.zeros([N,3],float)
    func[:,0]=x[:,0]
    func[:,1]=abs(phi[:,0])*2*(1e-9)
    return func

#Define the function that will calculate the relevant 'v' vector argument for the banded method. 
def v(phi,bv1,bv2):
    v=np.zeros([N+1,1],complex)
    v[0,0]=phi[0,0]*bv1[0,0]+bv2[1,0]*phi[1,0]
    for i in range(1,N):
        v[i,0]=bv1[i,0]*phi[i,0]+(bv2[i+1,0]*phi[i+1,0])+(bv2[i-1,0]*phi[i-1,0])                                     
    v[-1,0]=bv1[-1,0]*phi[-1,0]+bv2[N-1,0]*phi[N-1,0]
    return v

#Define function that will return array for the harmonic oscillator potential. 
def harmonicoscillator(x,k):
        L=len(x[:,0])
        V=[]
        for i in range(L):
                V.append(.5*k*(x[i,0]**2))
        return V

#Initialize the array for the 'A' tridiagonal matrix.
A=np.zeros([3,N+1],complex)

'''
Initializing arrays for the elements of the tridiagonal matricies.
Done to accomodate for the implemintatoin of varying potential
across the spatial domain.
'''

av1=np.zeros([N+1,1], complex)
av2=np.zeros([N+1,1], complex)
bv1=np.zeros([N+1,1], complex)
bv2=np.zeros([N+1,1], complex)

#Initialize x-position array
x=np.zeros([N+1,1],float)

#Creating spatial domain
for il in range(0,N+1):
    x[il,0]=(-L/2)+il*a


#Defining potential
V=harmonicoscillator(x,K)

'''
With potential defined and associated with respective x positions:
the elements of the tridiagonal matricies are constructed.
'''

for i in range(0,N+1):
    av1[i,0]=(h/2)*((2/h)+(1j*h_b/(m*(a**2)))+(1j*V[i]/h_b))
    av2[i,0]=-1j*h*h_b/4/m/(a**2)
    bv1[i,0]=-(h/2)*((1j*h_b/m/(a**2))+(1j*V[i]/h_b)-(2/h))
    bv2[i,0]=1j*h*h_b/4/m/(a**2)

#Filling the A tridiagonal matrix with coefficints, the structure is condusive to banded function. 
A[0,:],A[1,:],A[2,:]=av2[:,0],av1[:,0],av2[:,0]

def ClassicalPos(phi,x):
    EL=np.argmax(abs(phi[:,0]))
    pos=x[EL,0]
    return pos

#Position function for showing the potential curve in the region. 
def pos2(x,V):
    N=len(V)
    position=np.zeros([N,3],float)
    position[:,0]=x[:,0]
    position[:,1]=V
    position[:,1]*=1e7
    return position
phi=psi_i(xi,xf,x0,sig,k,N)

#Initialize phi_prime.
phi_prime=np.zeros([1,N+1],complex)

#Initialize x-position array
x=np.zeros([N+1,1],float)

#Fill x-position array. 
for il in range(0,N+1):
    x[il,0]=(-L/2)+il*a
    
#Create visual representation potential curve, initialize the curve object.

function=curve(pos=pos(phi,x),color=color.black)
Potential=curve(pos=pos2(x,V),color=color.blue)



t=0
while (True):
    rate(1e14)
    function.pos=pos(phi,x)
    phi=banded(A,v(phi,bv1,bv2),1,1)
    

