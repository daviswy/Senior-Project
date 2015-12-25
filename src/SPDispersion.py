#Wyatt Davis
#Final Project (M.N. Computational Physics, problem 9.8)

#Import all relevant functions.
import numpy as np  
import math as math
import cmath as cmath
from banded import banded
#import matplotlib as matplot
import matplotlib.pyplot as plt
import scipy.fftpack as fft

#Prompt user to choose a routine to run. 
E=input("Input a value for the energy (ev) for the incident particle: ")
p=input("Input 0 to test simulation for half period, 1 to run indefinitely, or 2 to see function after significant time evolution: ")

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


    
def find(max,phi):
    L=len(phi[:,0])
    for i in range(L-1):
        if phi[i,0]==max:
            return i

def dx(phi,x,max1,max2,a):
    half=max1/2.
    tol=.01
    L=len(phi[:,0])
    for i in range(L-1):
        if abs(phi[i,0])-half<tol:
            d1=x[i,0]
    dx=abs(d1-a*find(max1,abs(phi)))
    return dx

kx=[]
dk=1/N/a

for i in range(N//2+1):
    kx.append(i*dk)

def dft(y):
    N=len(y)
    c=np.zeros(N//2+1,complex)
    #loop that iterates through indicies representative of the frequency
    for k in range (N//2+1):
        '''loop that iterates through indicies representative of the 'y'
    data set values and calulates the relevant sum'''
        for n in range(N):
            c[k]+=y[n]*cmath.exp(-1j*2*math.pi*n*k/N)
    return c

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
#Initialize x-position array
x=np.zeros([N+1,1],float)

#Fill x-position array. 
for il in range(0,N+1):
    x[il,0]=(-L/2)+il*a

if p==0:
    
    while t<T/2:
        rate(1e14)
        phi=banded(A,v(phi,b1,b2),1,1)
        t+=h
        print(abs(np.max(phi)))
        
elif p==1:
    E*=1.602176e-19
    k=math.sqrt(2*m*E/h_b**2)
    phi=psi_i(xi,xf,x0,sig,k,N)
    M=[]
    t=[]
    xu=[]
    Max=[]
    Max2=.999858759542
    #Begin updating function position.
    time=0
    for i in range(100):
        phi=banded(A,v(phi,b1,b2),1,1)
        Max1=np.max(abs(phi))
        M.append(Max1)
        t.append(time)
        time+=1
        xu.append(dx(phi,x,Max1,Max2,a))
        Max2=Max1
        #s=dft(phi[:,0])

    



#print t 
#print M
plt.plot(t,xu)
plt.plot()
plt.show()



