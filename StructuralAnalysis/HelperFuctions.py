from math import inf, cos, pi, acos


def getQubicRoots(a,b,c,d):

    p = c-b**2/3
    q = d-(b*c)/3+(2*b**3)/27
    Delta = -4*p**3-27*q**2

    if Delta == 0:
        r1 = (-4*q)**(1/3)-b/3
        r2 = (q/2)**(1/3)-b/3
        r3 = r2
    elif Delta < 0:
        z1 = (-q/2+((p/3)**3+(q/2)**2)**0.5)**(1/3)
        z2 = (-q/2-((p/3)**3+(q/2)**2)**0.5)**(1/3)
        r1 = z1+z2-b/3
        r2 = inf
        r3 = inf
    else:
        thata = 1/3*acos(-q/2*(3/-p)**(3/2))
        r1 = 2*(-p/3)**0.5*cos(thata)-b/3
        r2 = 2*(-p/3)**0.5*cos(thata+2*pi/3)-b/3
        r3 = 2*(-p/3)**0.5*cos(thata+4*pi/3)-b/3

    return r1, r2, r3