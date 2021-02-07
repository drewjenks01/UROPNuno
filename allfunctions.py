# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 16:33:09 2021

@author: drewj
"""

import math

eps= 8.854*10**-12   #C^2/(N*m^2)
k= 1/(4*math.pi*eps)   #kg⋅m³⋅s⁻²⋅C^-2
mass= 9.109*10**-31     #kg
elecCharge=-1.602*10**-19  #C

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


    