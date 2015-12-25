#Wyatt Davis
#Senior Project

#Import all relevant functions.
import numpy as np  
import math as math
import cmath as cmath
from visual import curve,rate,color,arrow,text
from banded import banded


#Prompt user to choose a routine to run. 
p=input("Input 0 to test simulation for half period, 1 to run indefinitely, or 2 to see function after significant time evolution: ")
#prompt user to input energy values for incident particle and 
E=input("Input a value for the energy (ev) for the incident particle(E>50ev recommended for sake of speed): ")
EV=input("Input a value for the energy (ev) of the potential barrier: ")
EV*=1.602176e-19
E*=1.602176e-19
#Assign a variable to the length of the well. 
L=1e-8
#Assign a variable to the sigma value.
sig=1e-10
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
x0=L/2.5
#Assign a value to hbar. 
h_b=1.053e-34
#Assign a value to the mass of the electron.
m=9.109e-31
'''
Assign a variable to the k value for wave packet representation of particle,
utilizing classical relation between classical momentum and energy in addition
to wave number k relation to momentum.
'''
k=math.sqrt(2*m*E/h_b**2)

#assign value to constant potential
V0=EV
k1=math.sqrt(2*m*(E-V0)/h_b**2)

# This function calculates the t=0 array for psi for all x between 0 and L. Keeps x=0,L at zero in order to satisfy the boundary conditions. 
def psi_i(xi,xf,x0,sig,k,N):
	a=(xf-xi)/N
	psi=np.zeros([N+1,1], complex)
	for i in range(0,N):
		psi[i,0]=cmath.exp(1j*k*i*a)
	return psi

#Define the position function that I will utilize in updating the position of the curve object. 
def pos(phi,x):
    N=len(phi)
    func=np.zeros([N,3],float)
    func[:,0]=x[:,0]
    func[:,1]=abs(phi[:,0])*4*(1e-9)
    return func

#Define the function that will calculate the relevant 'v' vector argument for the banded method. 
def v(phi,b1,b2,bv1,bv2):
    v=np.zeros([N+1,1],complex)
    v[0,0]=phi[0,0]*b1+b2*phi[1,0]
    for i in range(1,(N/2)):
        v[i,0]=b1*phi[i,0]+b2*(phi[i+1,0]+phi[i-1,0])
    
    for n in range(((N/2)),N-1):
        v[n,0]=bv1*phi[n,0]+bv2*(phi[n+1,0]+phi[n-1,0])
                                         
    v[-1,0]=bv1*phi[-1,0]+bv2*phi[N-1,0]
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



#Define values for region of constant potential
av1=(h/2)*((2/h)+(1j*h_b/(m*(a**2)))+(1j*V0/h_b))
av2=-1j*h*h_b/4/m/(a**2)
bv1=-(h/2)*((1j*h_b/m/(a**2))+(1j*V0/h_b)-(2/h))
bv2=1j*h*h_b/4/m/(a**2)

#Filling the A-array with constants, structure condusive to banded function. 
A[0,0:N/2],A[1,0:N/2],A[2,0:N/2]=a2,a1,a2
A[0,N/2:(N+1)],A[1,N/2:(N+1)],A[2,N/2:(N+1)]=av2,av1,av2    
#Initialize phi at t=0.

phi=psi_i(xi,xf,x0,sig,k,N)

#Initialize phi_prime.
phi_prime=np.zeros([1,N+1],complex)

#Initialize x-position array
x=np.zeros([N+1,1],float)

#Fill x-position array. 
for il in range(0,N+1):
    x[il,0]=(-L/2)+il*a
kx=[]
dkx=1/a/N
for g in range(N//2+1):
    kx.append(g*dkx)

R=((k-k1)/(k+k1))**2
print(R)
if p==0:
    '''
    Create visual representation of box potential, initialize the curve object,
    and a reference arrow at the initial position
    '''
    center=arrow(pos=(0,0,0),axis=(0,4e-9,0),shaftwidth=5e-11,color=color.white)
    function=curve(pos=pos(phi,x),color=color.red)
    Wall1=arrow(pos=(-L/2,0,0),axis=(0,4e-9,0),shaftwidth=10e-11, color=color.blue)
    Wall2=arrow(pos=(L/2,0,0),axis=(0,4e-9,0),shaftwidth=10e-11,color=color.blue)
    #Initialize time. 
    t=0
    #Begin updating function position, this loop will break after it reaches a half cycle. 
    while t<T/2:
        rate(1e14)
        function.pos=pos(phi,x)
        phi=banded(A,v(phi,b1,b2),1,1)
        t+=h
        print(abs(np.max(phi)))




        
        
elif p==1:
    #Create visual representation of box potential, initialize the curve object.
    function=curve(pos=pos(phi,x),color=color.red)
    Potential=curve(pos=[(-L/2,-.2e-9,0),(0,-.2e-9,0),(0,(EV/(E+(EV)))*5e-9,0),(L/2,(EV/(E+EV))*5e-9,0)], color=color.blue)
    EnLab=text(text=repr(E/1.602176e-19) , align='center',pos=(-L/4,((E/(E+EV))*5e-9)+.5e-9,0),height=.25e-9,depth=.1e-9,color=color.white)
    PoLab=text(text=repr(EV/1.602176e-19), align='center',pos=(L/4,((EV/(E+EV))*5e-9)+.5e-9,0),height=.25e-9,depth=.1e-9,color=color.blue)
    EnergyI=curve(pos=[(-L/2,(E/(E+EV))*5e-9,0),(L/2,(E/(EV+E))*5e-9,0)],color=color.white)
    Wall1=arrow(pos=(-L/2,0,0),axis=(0,4e-9,0),shaftwidth=10e-11, color=color.blue)
    Wall2=arrow(pos=(L/2,0,0),axis=(0,4e-9,0),shaftwidth=10e-11,color=color.blue)
    #Begin updating function position.
    t=0
    #Normalizer=(np.trapz(np.absolute(phi[0:N/2,0]),dx=a))**2
    
    xf=np.linspace
    while True:
        rate(1e14)
        function.pos=pos(phi,x)
        phi=banded(A,v(phi,b1,b2,bv1,bv2),1,1)
        t+=h
    
        #print((np.trapz(np.absolute(phi[0:N/2,0]),dx=a)**2)/(np.trapz(np.absolute(phi[:,0]),dx=a)**2),R)
        #print(np.trapz(np.absolute(phi[0:N/2,0]),dx=a))
        print((np.trapz(np.absolute(phi[0:N/2,0])**2,dx=a))/(np.trapz(np.absolute(phi[:,0])**2,dx=a)),R)

elif p==2:
    #Prompt user to input number of periods into animation they would like to start animation. 
    P=input("Enter the number of classical periods into time evolution you would like to see the animation start: ")
    Time=T*P
    print Time
    t=0
    #This loop updates phi but without the added animation. Speeding up the calculation slightly.
    while Time>t:
        phi=banded(A,v(phi,b1,b2),1,1)
        t+=h
    #Create visual representation of box potential, initialize the curve object.
    Wall1=arrow(pos=(-L/2,0,0),axis=(0,4e-9,0),shaftwidth=10e-11, color=color.blue)
    Wall2=arrow(pos=(L/2,0,0),axis=(0,4e-9,0),shaftwidth=10e-11,color=color.blue)
    function=curve(pos=pos(phi,x),color=color.red)
    #Begin updating function position.
    while True:
        rate(1e14)
        function.pos=pos(phi,x)
        phi=banded(A,v(phi,b1,b2),1,1)
        


