# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 11:36:32 2021

@author: drewj
"""

import math
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from random import uniform
import random
from timeit import default_timer as timer
from scipy import integrate



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
        return [self.getPos(),self.getVel()]
    
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
        dvdt= electron.getFor()/mass
        dxdt=z[1]
        return[dxdt,dvdt]

def elecPot(start,end,efield):
    return efield*(end-start)

def elecPot2(pos,efield):
    return Efieldleft*Eleftbot-efield*pos+Efieldright*Erightbot
    

########################################
########   INITIAL SETTINGS    #########
########################################

#sets the timepoints
n=60000
t=np.linspace(0,2,n)

#number of electrons in system
num_e = 2

#list for initial positions
initpos=[-2,2]  

#creats random starting positions for each electrons b/w a certain range
# while(len(initpos)<num_e):
#     for i in range(num_e):
#               r=uniform(-30,30)
#               if r not in initpos: initpos.append(r)
# print(initpos)

#list of electrons
electrons = [] 

Efieldleft=-10**-7 #V/m
Efieldright=10**-7
Elefttop=-10
Eleftbot=-8
Erighttop=10
Erightbot=8

 #creates num_e number of electrons w/ rand pos and 0 vel
for i in range (num_e): 
    x0 = initpos[i]
    v0 = 0
    f0 = 0
    electrons.append(electron(x0,v0,f0)) #adds new instance of electron to an array
    #print(x0)
    
    
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
                
        
        if(electron.getPos()>Erightbot):
             fnew= Efieldright*elecCharge
             pe_now-=elecPot2(electron.getPos()-Erightbot,Efieldright)*elecCharge
        if(electron.getPos()<Eleftbot):
             fnew= Efieldleft*elecCharge
             pe_now+=elecPot2(electron.getPos()-Eleftbot,Efieldleft)*elecCharge
               
                
        fnet+=fnew
            
       # print("2",fnet)
                
        electron.ChangeFor(fnet)
    
    pe[x] = pe_now
    totale[x] = pe[x] + ke[x]

    for electron in electrons:
       # print(electron.getVel())
        z = odeint(model1,electron.getPosVel(),t)
        electron.ChangePos(z[1][0],x)
        electron.ChangeVel(z[1][1],x)
       # print(electron.getVel())
        
    #print(totale[x])
    
    #print([x,electron.getPosListi(x),electron.getVelListi(x)])
        

end=timer()
print(end-start)
#print(electrons[1].getPos())


#########################################
########       PLOTTING         #########
#########################################

#time scale data
scale=[6,9,12,16,19,24,26.5,30.5,36,75,280,753,1790,3383]
num=[2,3,4,5,6,7,8,9,10,20,50,100,200,300]
data = [] #example
for i in range(len(num)):
    data.append(num[i]**2)

#time scale plot
# plt.plot(num,scale)
# plt.plot(num,data)
# plt.xlabel('N')
# plt.ylabel('time')
# plt.yscale('log')
# plt.xscale('log')
# #plt.axis([0,2,-0.5,0.5])
# plt.legend(['real data','example'])
# plt.show()
    
#position vs time
for electron in electrons:
    plt.plot(t,electron.getPosList())
    plt.xlabel('time')
    plt.ylabel('position')
plt.axis([1.5,2,-50,50])
plt.show()

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