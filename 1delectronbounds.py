# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 16:32:28 2021

@author: drewj
"""
import electronclass as ec
import allfunctions as af
import math
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from random import uniform
from timeit import default_timer as timer



#########################################
########    CONSTANTS           #########
#########################################

eps= 8.854*10**-12   #C^2/(N*m^2)
k= 1/(4*math.pi*eps)   #kg⋅m³⋅s⁻²⋅C^-2
mass= 9.109*10**-31     #kg
elecCharge=-1.602*10**-19  #C


#########################################
########    ODE Model           #########
#########################################
def model1(z,t):
        return[z[1],electron.getFor()/mass]
    
    
########################################
########   INITIAL SETTINGS    #########
########################################

#sets the timepoints
n=40000
t=np.linspace(0,4,n)

#number of electrons in system
num_e = 3

#list for initial positions
initpos=[]  

#lists of electrons
electrons,staticelectrons = [],[]

#creats random starting positions for each electrons b/w a certain range
while(len(initpos)<num_e):
    for i in range(num_e):
              r=uniform(-10,10)
              if r not in initpos: 
                  initpos.append(r)
                  electrons.append(ec.electron(r,0,0,n)) #adds new instance of electron to an array


staticelectrons.append(ec.electron(20,0,0,n))
staticelectrons.append(ec.electron(-20,0,0,n))
    
    
#creates arrays for energy values
ke, pe, totale = np.zeros(n),np.zeros(n),np.zeros(n)


#starts timer to check how long the main loop takes
start=timer()


########################################
########      MAIN LOOP        #########
########################################

for x in range (0,n):
    ke_now,pe_now,counter = 0,0,0
    pos_now = []
    
    for electron in electrons:      #makes an array of current electron positions and calculates KE
        pos_now.append(electron.getPos())
        ke_now +=af.Kinetic(electron.getVel())  
    
    for electron in electrons:      #solves for PE
        for i in range (counter+1,num_e):
            pedist = af.getDistance(electron.getPos(), pos_now[i])
            pe_now +=af.PotE(abs(pedist))
        counter +=1  
        pe_now+=af.PotE(abs(af.getDistance(electron.getPos(),staticelectrons[0].getPos())))+af.PotE(abs(af.getDistance(electron.getPos(),staticelectrons[1].getPos())))

        fnet = 0
        for i in range (num_e):     #solves for total force on each electron
            distance = af.getDistance(electron.getPos(),pos_now[i])
            if (distance==0):
                fnet += 0
            elif (distance>0):
                fnet += k*(elecCharge**2)/(distance**2) #force will push electron right
            else:
                fnet -= k*(elecCharge**2)/(distance**2) #force will push electron left
                
        fnet+= -k*(elecCharge**2)/(af.getDistance(electron.getPos(),staticelectrons[0].getPos())**2)+k*(elecCharge**2)/(af.getDistance(electron.getPos(),staticelectrons[1].getPos())**2)
        electron.ChangeFor(fnet)
        z = odeint(model1,electron.getPosVel(),t)   #gets the new pos and vel
        electron.ChangePos(z[1][0],x)       #updates pos
        electron.ChangeVel(z[1][1],x)       #updates vel
        
    ke[x],pe[x] = ke_now,pe_now         #adds the new KE and PE to arrays
    totale[x]=pe[x] + ke[x]             #addds total energy to array
    
end=timer()
print(end-start)



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