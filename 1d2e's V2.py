# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 14:41:08 2021

@author: drewj
"""

import math
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
from random import randrange

#########################################
########    CONSTANTS           #########
#########################################

NumElec= 2   #unused so far
eps= 8.854*10**-12   #C^2/(N*m^2)
k= 1/(4*math.pi*eps)   #kg⋅m³⋅s⁻²⋅C^-2
mass= 9.109*10**-31     #kg
elecCharge=-1.602*10**-19  #C

#########################################
########    CREATE ELECTRONS    #########
#########################################

class electron():
    
    def __init__(self,x,v):
        self.x=x
        self.v=v
        self.poslist=np.zeros(n)
        self.vellist=np.zeros(n)
       
    def ChangePos(self,xnew,i):     #sets the new current position of the electron
        self.x=xnew
        self.AddPos(xnew,i)
    
    def ChangeVel(self,vnew,i):     #sets the new current vel of the electron
        self.v=vnew
        self.AddVel(vnew,i)
    
    def getPos(self):       #return the electron's position
        return self.x
    
    def getVel(self):       #returns the electron's velocity
        return self.v
    
    def AddPos(self,value,i):       #adds positions to the position array
       self.poslist[i]=value
       
    def AddVel(self,value,i):       #adds velocities to the velocity array
       self.vellist[i]=value
       
    def getPosListi(self,i):        #returns a sepecific index from the position array
        return self.poslist[i]
    
    def getVelListi(self,i):         #returns a sepecific index from the velocity array
        return self.vellist[i]
    
    def getPosList(self):       #returns the position array
        return self.poslist
    
    def getVelList(self):       #returns the velocity array
        return self.vellist
    
    def getPosVel(self):        #returns the current position and velocity as an array
        return [self.getPos(),self.getVel()]
    
    def CheckPos(self,pos2):        #checks if the electrons is to the left or right of another position
        if(self.x-pos2>0):
            return 1
        else:
            return 0
        
        
#########################################
########    FUNCTIONS           #########
#########################################

#finds the distance from some other position
def getDistance(d1,d2):       
        return abs(d1-d2)
    
#solves ODE with force due to boundary distance as negative
#solves ODE with force due to boundary distance as negative
def model1(z,t):
    x1=d1  #sets x1 equal to distance b/w middle
 #   v=z[1] #sets v equal to 2nd array position
    force=k*(elecCharge**2)/((x1**2)*mass)
    dvdt= force
    dxdt=z[1]
    return[dxdt,dvdt]

#solves ODE with force due to distance b/w middle as negative
def model2(z,t):
    x1=d1
    v=z[1]
    force=-k*(elecCharge**2)/((x1**2)*mass)
    dvdt= force
    dxdt=v
    return[dxdt,dvdt]

#calculates potential energy
def ElecPot(k,Q,dist):
    return (k*Q**2)/(dist)

#calculates kinetic energy
def Kinetic(mass,vel):
    return (1/2)*mass*(vel**2)
        
        

#sets the timepoints
n=10000
t=np.linspace(0,5,n)

#creates electrons

e1=electron(1,0)        #starts at x=1 and v=0
e2=electron(-1,0)       #starts at x=1 and v=0

#adds initial values to respective arrays
e1.AddPos(e1.getPos(),0)
e2.AddPos(e2.getPos(),0)

#create arrays to hold values
ke=np.zeros(n)  #kinetic energy
pe=np.zeros(n)  #potential energy
energy=np.zeros(n)  #total energy


#inital energy values
ke[0]=0     
pe[0]=ElecPot(k,elecCharge,2)
energy[0]=ke[0]+pe[0]


#solves the ODEs, records the outputs, and moves forward to the next timestep, repeat
for i in range(1,n):
    
    #find the distance
    d1=e1.getPos()-e2.getPos()
   
    #solve ODE's
    z= odeint(model1,e1.getPosVel(),t)
   # print(z)
    r= odeint(model2,e2.getPosVel(),t)

    #set new positions and velocities
    z0=z[1]
    e1.ChangePos(z0[0],i)
    e1.ChangeVel(z0[1],i)
    
    r0=r[1]
    e2.ChangePos(r0[0],i)
    e2.ChangeVel(r0[1],i)
    
    #finds the new distance (only useful for the energy calculations)
    d1new=e1.getPos()-e2.getPos()
    
    #find KE,PE,TE and add to arrays
    ke[i]=Kinetic(mass,e1.getVel())+Kinetic(mass,e2.getVel())
    pe[i]=abs(ElecPot(k,elecCharge,d1new))
    energy[i]=ke[i]+pe[i]
  #  print(energy[i]-energy[i-1])
    
        
  #  print timesteps, positions, and velocities for both electrons
    if i==1:  
        print(["Timestep","Position","Velocity"])
        print([0,e1.getPosListi(0),e1.getVelListi(0)])
        print([0,e2.getPosListi(0),e2.getVelListi(0)])
    print([i,e1.getPosListi(i),e1.getVelListi(i)])
    print([i,e2.getPosListi(i),e2.getVelListi(i)])
    


#########################################
########       PLOTTING         #########
#########################################

#plot positions
plt.plot(t,e1.getPosList(),'r-')
plt.plot(t,e2.getPosList(),'g-')
plt.xlabel('time')
plt.ylabel('position')
plt.legend(['e1','e2'])
plt.show()

#plot velocities
# plt.plot(t,e1.getVelList(),'r--')
# plt.plot(t,e2.getVelList(),'g--')
# plt.xlabel('time')
# plt.ylabel('velocity')
# plt.legend(['e1','e2'])
# plt.show()

# #plot energies
# plt.plot(t,energy,'g-')
# plt.plot(t,pe,'r-')
# plt.plot(t,ke,'b-')
# plt.xlabel('time')
# plt.ylabel('energy')
# plt.legend(['total','PE','KE'])
# #plt.axis([-0.5,60,0,(1.5)*(10**-28)])
# plt.show()


# #create velocities from positions
# #vels[i]=math.sqrt(((k*elecCharge**2)/mass)*((1/2)-(1/getDistance(e1.getPosListi(i),e2.getPosListi(i)))))
  
    

# #plot velocity vs position
# plt.plot(e1.getPosList(),e1.getVelList(),'r-')   #based off of our code
# plt.plot(e1.getPosList(),vels,'b-')
# #plt.plot(e1.getPosList(),vels2,'g-',alpha=0.7)
# plt.xlabel('position')
# plt.ylabel('velocity')
# plt.legend(['code','equation','equation2'])
# plt.show()

    
    
    
    
    
    
    
    
    
    
    
    
    
    