

#Wyatt Davis
#Final Project (M.N. Computational Physics, problem 9.8)

#Import all relevant functions.
import numpy as np  
import math as math
import cmath as cmath
from visual import curve,rate,color,arrow,display 
from banded import banded

scene=display(title='Wave Packet Evolution: Free Particle', background=color.white)

#Prompt user to choose a routine to run. 
p=input("Input 0 to test simulation for half period, 1 to run indefinitely, or 2 to see function after significant time evolution: ")

#Assign a variable to the length of the well. 
L=1e-8
#Assign a variable to the sigma value.
sig=1e-10
#Assign a variable to the k value. 
k=0
#Assign a variable to the initial x position.
xi=0.
#Assign a variable to the length of the well.
xf=L
#Assign a variable to number of intervals.
N=1000
#Assign a variable to the grid spacing.
a=(xf-xi)/N
#Assign a variable to the time step.
h=1e-18
#Assign a value to the initial x position.
x0=L/5
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
    func[:,1]=phi[:,0]*3*(1e-9)
    return func

#Define the function that will calculate the relevant 'v' vector argument for the banded method. 
def v(phi,b1,b2):
    v=np.zeros([N+1,1],complex)
    v[0,0]=phi[0,0]*b1+b2*phi[1,0]
    for i in range(1,N-1):
        v[i,0]=b1*phi[i,0]+b2*(phi[i+1,0]+phi[i-1,0])
    v[-1,0]=b1*phi[-1,0]
    return v 

#Initialize the array.
A=np.zeros([3,N+1],complex)
'''
Define the values to go into the banded matricies: 'A' and 'v'.
These are constants pulled from the update equation
'''
a1=1+h*(1j*h_b/2/m/(a**2))
a2=-h*(1j*h_b/4/m/(a**2))
b1=1-h*(1j*h_b/2/m/(a**2))
b2=h*(1j*h_b/4/m/(a**2))

#Filling the A-array with constants, structure condusive to banded function. 
A[0,:],A[1,:],A[2,:]=a2,a1,a2
    
#Initialize phi at t=0. 
phi=psi_i(xi,xf,x0,sig,k,N)

#Initialize phi_prime.
phi_prime=np.zeros([1,N+1],complex)

#Initialize x-position array
x=np.zeros([N+1,1],float)

#Fill x-position array. 
for il in range(0,N+1):
    x[il,0]=(-L/2)+il*a





'''
The following loop will calculate the new phi value via
the banded method (written by Mark Newman), which updates the position of the curve
object: 'function', by solving for all of the new array points at the same time
though a set of equations in matrix form.I update the curve object position by
using a user defined function: 'pos', which
takes the phi and x-position arrays as arguments and returns a relevant 2D
array. After each iteration the 't' phi array is replaced with the 't+h' phi array
allowing for the updating of the curve object with time.  
'''

        


function=curve(pos=pos(phi,x),color=color.black)
    
while True:
        rate(1e14)
        function.pos=pos(phi,x)
        phi=banded(A,v(phi,b1,b2),1,1)
        
        p=np.argmax(phi[:,0])
        print(abs(phi[p,0]))


