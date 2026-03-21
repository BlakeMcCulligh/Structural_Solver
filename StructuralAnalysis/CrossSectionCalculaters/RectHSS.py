
def getA(d, b, t):
    return b*d-(b-2*t)*(d-2*t)

def getIx(d, b, t):
    return 1/12*(b*d**3-(b-2*t)*(d-2*t)**3)

def getIy(d, b, t):
    return 1/12*(d*b**3-(d-2*t)*(b-2*t)**3)

def getJ(d, b, t):
    return 2*t*d**2*b**2/(d+b)