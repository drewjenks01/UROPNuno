# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 17:34:16 2021

@author: drewj
"""

import math
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random
from random import randrange
#import matplotlib.animation as animation

#constants
NumElec= 2   #unused so far
eps= 8.854*10**-12   #C^2/(N*m^2)
k= 1/(4*math.pi*eps)   #kg⋅m³⋅s⁻²⋅C^-2
mass= 9.109*10**-31     #kg
elecCharge=-1.602*10**-19  #C

#create electrons
class electron():
    
    def __init__(self,x,v,f):
        self.x=x
        self.v=v
        self.f=f
        self.poslist=np.zeros(n)
        self.vellist=np.zeros(n)

       
    def ChangePos(self,xnew,i):     #sets the new current position of the electron
        self.x=xnew
        self.AddPos(xnew,i)
    
    def ChangeVel(self,vnew,i):     #sets the new current vel of the electron
        self.v=vnew
        self.AddVel(vnew,i)
    
    def ChangeFor(self,fnew):
        self.f=fnew
 
    def getFor(self):
        return self.f
    
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
        

def getDistance(d1,d2):       #finds the distance from some other position
    return (d1-d2)

#calculates potential energy
def ElecPot(k,Q,dist):
    return (k*Q**2)/dist

#calculates kinetic energy
def Kinetic(mass,vel):
    return (1/2)*mass*(vel**2)

#sets the timepoints
n=10000
t=np.linspace(0,5,n)

num_e = 3 #number of electrons in system

initpos=[]

for i in range(0,num_e):
    rand=randrange(-10,10)    #makes a random int from -10-10
    initpos.append(rand)
    
electrons = [] #list of electrons

for i in range (num_e):  #creates num_e number of electrons w/ rand pos and 0 vel
    x0 = initpos[i]
    v0 = 0
    f0 = 0
    electrons.append(electron(x0,v0,f0)) #adds new instance of electron to an array
    #print(x0)

def model1(z,t):
        v=z[1] #sets v equal to 2nd imputed array position
        dvdt= electron.getFor()/mass
        dxdt=v
        return[dxdt,dvdt]
    
ke = np.zeros(n)
pe = np.zeros(n)
totale = np.zeros(n)

ke[0] = 0
pe[0] = 0
totale[0] = ke[0] + pe[0]

for x in range (0,n):
    ke_now = 0
    pos_now = []
    
    for electron in electrons:
        pos_now.append(electron.getPos())
        ke_now = ke_now+Kinetic(mass, electron.getVel())  
        
    ke[x] = ke_now
    
   
    pe_now = 0
    counter = 0
    for electron in electrons:
        for i in range (counter,num_e):
            pedist = getDistance(electron.getPos(), pos_now[i])
            if (pedist != 0):
                pe_now = pe_now + ElecPot(k, elecCharge, abs(pedist))
        counter = counter +1           
    
    pe[x] = pe_now
    totale[x] = pe[x] + ke[x]
    
    for electron in electrons:
        fnet = 0
        distances = np.zeros(num_e)
        for i in range (num_e):
            distances[i] = getDistance(electron.getPos(),pos_now[i])
        #print(distances)
        
        for i in range (num_e):
            if (distances[i]==0):
                force = 0
            if (distances[i]>0):
                force = k*(elecCharge**2)/((distances[i]**2)) #force will push electron right
            if (distances[i]<0):
                force = -k*(elecCharge**2)/((distances[i]**2)) #force will push electron left
            
            fnet = fnet + force
                
        electron.ChangeFor(fnet)
        #fnets.append(fnet)
        #print(fnet)
        #print(electron.getFor())
    
    for electron in electrons:
        z = odeint(model1,electron.getPosVel(),t)
        z0=z[1]
        electron.ChangePos(z0[0],x)
        electron.ChangeVel(z0[1],x)
        
    #print(totale[x])
    
    #print([x,electron.getPosListi(x),electron.getVelListi(x)])
        
        
 
for electron in electrons:
    plt.plot(t,electron.getPosList())
    plt.xlabel('time')
    plt.ylabel('position')

plt.show()
    
for electron in electrons:
    plt.plot(t,electron.getVelList())
    plt.xlabel('time')
    plt.ylabel('velocity')

plt.show()


plt.plot(t,totale,'g-')
plt.plot(t,pe,'r-')
plt.plot(t,ke,'b-')
plt.xlabel('time')
plt.ylabel('energy')
plt.legend(['total','PE','KE'])
#plt.axis([-0.5,60,0,(1.5)*(10**-28)])
plt.show()