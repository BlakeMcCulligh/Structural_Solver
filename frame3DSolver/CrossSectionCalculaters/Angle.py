
def getX(d, b, t):
    return (b**2+(d-t)*t)/(2*(b+d-t))
def getY(d, b, t):
    return (d**2+(b-t)*t)/(2*(b+d-t))

def getA(d, b, t):
    return t*(b+d-t)

def getIx(d, b, t):
    y = getY(d, b, t)
    return (t*(d-y)**3+b*y**3-(b-t)*(y-t)**3)/3
def getIy(d, b, t):
    x = getX(d, b, t)
    return (t*(b-x)**3+d*x**3-(d-t)*(x-t)**3)/3

def getJ(d, b, t):
    # assuming t is small compared to d and b
    return t**3*(b+d)/3