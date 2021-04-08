# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 22:58:05 2021

@author: drewj
"""

import math
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from random import uniform
from timeit import default_timer as timer



#import matplotlib.animation as animation

#########################################
########    CONSTANTS           #########
#########################################

eps= 8.854*10**-12   #C^2/(N*m^2)
k= 1/(4*math.pi*eps)   #kg⋅m³⋅s⁻²⋅C^-2
mass= 9.109*10**-31     #kg
elecCharge=-1.602*10**-19  #C

#########################################
########    CREATE ELECTRONS    #########
#########################################

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
        return [self.x,self.v]
    
    def getForVel(self):
        return [self.f/9.109*10**-31,self.v]
    
    def CheckPos(self,pos2):        #checks if the electrons is to the left or right of another position
        if(self.x-pos2>0):
            return 1
        else:
            return 0
        
#########################################
########       FUNCTIONS        #########
#########################################

def getDistance(d1,d2):       #finds the distance from some other position
    return (d1-d2)

#calculates potential energy
def PotE(dist):
    return (k*elecCharge**2)/dist

#calculates kinetic energy
def Kinetic(vel):
    return (1/2)*mass*(vel**2)

def model1(z,t):
        dvdt=electron.getFor()/mass
        dxdt=z[1]
        return[dxdt,dvdt]
    

########################################
########   INITIAL SETTINGS    #########
########################################

#sets the timepoints
n=20000
t=np.linspace(0,10,n)

#number of electrons in system
num_e = 3

#list for initial positions
initpos=[]  

#list of electrons
electrons = [] 
staticelectrons=[]  #boundary electrons

#creats random starting positions for each electrons b/w a certain range
while(len(electrons)<num_e):
    for i in range(num_e):
              r=uniform(-15,15)
              if r not in initpos: 
                  initpos.append(r)
                  electrons.append(electron(r,0,0))

#creates the static electrons acting as boundaries
staticelectrons.append(electron(20,0,0))
staticelectrons.append(electron(-20,0,0))
    
    
#creates arrays for energy values
ke = np.zeros(n)
pe = np.zeros(n)
totale = np.zeros(n)

#starts timer to check how long the main loop takes
start=timer()


########################################
########      MAIN LOOP        #########
########################################

for x in range (0,n):
    ke_now = 0
    pos_now = []
    
    for electron in electrons:
        pos_now.append(electron.getPos())
        ke_now +=Kinetic(electron.getVel())  
        
    ke[x] = ke_now
    
   
    pe_now = 0
    counter = 0
    for electron in electrons:
        for i in range (counter,num_e):
            pedist = getDistance(electron.getPos(), pos_now[i])
            if (pedist != 0):
                pe_now +=PotE(abs(pedist))
        counter +=1    
        pe_now+=PotE(abs(getDistance(electron.getPos(),staticelectrons[0].getPos())))
        pe_now+=PotE(abs(getDistance(electron.getPos(),staticelectrons[1].getPos())))
    pe[x] = pe_now
    totale[x] = pe[x] + ke[x]
    
    for electron in electrons:
        fnet = 0
        fnew=0
        distances = np.zeros(num_e)
        for i in range (num_e):
            distances[i] = getDistance(electron.getPos(),pos_now[i])
        
        for i in range (num_e):
            if (distances[i]==0):
                fnet += 0
            elif (distances[i]>0):
                fnet += k*(elecCharge**2)/((distances[i]**2)) #force will push electron right
            else:
                fnet += -k*(elecCharge**2)/((distances[i]**2)) #force will push electron left
                
        fnet-= k*(elecCharge**2)/(getDistance(electron.getPos(),staticelectrons[0].getPos())**2)
        fnet+= k*(elecCharge**2)/(getDistance(electron.getPos(),staticelectrons[1].getPos())**2)
        electron.ChangeFor(fnet)

        z = odeint(model1,electron.getPosVel(),t)
        electron.ChangePos(z[1][0],x)
        electron.ChangeVel(z[1][1],x)
        

end=timer()
print(end-start)
#print(electrons[1].getPos())


#########################################
########       PLOTTING         #########
#########################################


#position vs time
for electron in electrons:
    plt.plot(t,electron.getPosList())
    plt.xlabel('time')
    plt.ylabel('position')
# plt.axis([14,25,-55,55])
plt.show()
    
# #velocity vs time
for electron in electrons:
    plt.plot(t,electron.getVelList())
    plt.xlabel('time')
    plt.ylabel('velocity')
plt.show()


# #energy vs time
plt.plot(t,totale,'g-')
plt.plot(t,pe,'r-')
plt.plot(t,ke,'b-')
plt.xlabel('time')
plt.ylabel('energy')
plt.legend(['total','PE','KE'])
# #plt.axis([-0.5,60,0,(1.5)*(10**-28)])
plt.show()