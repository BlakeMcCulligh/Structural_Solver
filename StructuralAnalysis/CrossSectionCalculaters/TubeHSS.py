from math import pi

def getA(d, t):
    return pi*(d**2-(d-2*t)**2)/4

def getI(d, t):
    return pi*(d**4-(d-2*t)**4)/64

def getJ(d, t):
    return pi*(d**4-(d-2*t)**4)/32