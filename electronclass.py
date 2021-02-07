# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 16:33:33 2021

@author: drewj
"""
import numpy as np
#########################################
########    CREATE ELECTRONS    #########
#########################################

class electron():
    
    def __init__(self,x,v,f,n):
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
    
    def CheckPos(self,pos2):        #checks if the electrons is to the left or right of another position
        if(self.x-pos2>0):
            return 1
        else:
            return 0
    def getForVel(self):
        return [self.f/9.109*10**-31,self.v]