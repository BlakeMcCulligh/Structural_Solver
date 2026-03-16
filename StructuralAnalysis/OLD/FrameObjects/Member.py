from math import pi, acos, cos, inf

import numpy as np


class Member:
    def __init__(self):
        self.nodes = []
        self.restraints = [] # [[0,0,0],[0,0,0]]   0 = free, 1 = pin, 2 = moment
        self.crossSection = None
        self.material = None
        self.length = None
        self.K = None

    def calcLength(self):
        length = 0
        for i in range(self.nodes[0].numDimensions):
            length += (self.nodes[1].cords[i] - self.nodes[0].cords[i]) ** 2

        self.length = length ** 0.5

    def calcK(self):
        node1 = 2
        node2 = 2
        for i in range(self.nodes[0].numDimensions):
            if self.restraints[0][i] == 0:
                node1 = 0
            if self.restraints[1][i] == 0:
                node2 = 0

        for i in range(self.nodes[0].numDimensions):
            if self.restraints[0][i] == 1:
                node1 = 1
            if self.restraints[1][i] == 0:
                node2 = 1

        if node1 == 1 and node2 == 1:
            self.K = 1
        elif (node1 == 1 and node2 == 2) or (node1 == 2 and node2 == 1):
            self.K = 0.7
        elif node1 == 2 and node2 == 2:
            self.K = 0.5
        elif (node1 == 0 and node2 == 2) or (node1 == 2 and node2 == 0):
            self.K = 2


    def getTensionResistance(self):

        Tr = 0.9 * self.crossSection.Ag * self.material.Fy
        return Tr

    def getCompressionResistance(self):
        Fe = self.getFe()
        Lambda = (self.material.Fy / Fe) ** 0.5
        Cr = 0.9 * self.crossSection.A * self.material.Fy / (1 + Lambda ** (2 * self.material.n)) ** (1 / self.material.n)
        return Cr

    def getFe(self):
        ro2 = self.crossSection.xo ** 2 + self.crossSection.yo ** 2 + self.crossSection.rx ** 2 + self.crossSection.ry ** 2

        Fex = (pi ** 2) * self.material.E / (self.K * self.length / self.crossSection.rx)
        Fey = (pi ** 2) * self.material.E / (self.K * self.length / self.crossSection.ry)
        Fez = ((pi ** 2 * self.material.E * self.crossSection.Cw) / (self.K * self.length) ** 2 + self.material.G * self.crossSection.J) * 1 / self.crossSection.A * ro2

        if self.crossSection.numSymmetry == 2:
            return min(Fex, Fey, Fez)
        elif self.crossSection.numSymmetry == 1:
            ohm = 1 - (self.crossSection.xo ** 2 + self.crossSection.yo ** 2 / ro2)
            Feyz = (Fey + Fez) / (2 * ohm) * (1 - (1 - (4 * Fey * Fez * ohm) / (Fey + Fez) ** 2) ** 0.5)
            return min(Fex, Feyz)
        else:
            a = 1 - self.crossSection.xo ** 2 / ro2 - self.crossSection.yo ** 2 / ro2
            b = Fey * (self.crossSection.xo ** 2 / ro2) + Fex * (self.crossSection.yo ** 2 / ro2) - Fex - Fey - Fez
            c = Fey * Fez + Fex * Fez + Fex * Fey
            d = -1 * Fex * Fey * Fez

            [r1, r2, r3] = getQubicRoots(a, b, c, d)
            Fe = min(r1, r2, r3)
            return Fe

    def localStifnessMatrix(self):
        a=1
        # TODO
        # if self.crossSection.numSymmetry == 2:
        #     a = 1
        #
        #
        # elif self.nodes[0].numDimensions == 3:
        #     # 3D local stiffness matix including shear and bending
        #     Constant_y = (12 * self.material.E * self.crossSection.Iz) / (
        #                 self.material.G * self.material.Asy * self.length ** 2)
        #     Constant_z = (12 * self.material.E * self.crossSection.Iy) / (
        #             self.material.G * self.material.Asz * self.length ** 2)
        #
        #     T = self.material.E * self.crossSection.Ax / self.length
        #
        #     bend12Z = (12 * self.material.E * self.crossSection.Iz) / (self.length ** 3 * (1 + Constant_y))
        #     bend12Y = (12 * self.material.E * self.crossSection.Iy) / (self.length ** 3 * (1 + Constant_z))
        #     twisting = (self.material.G*self.crossSection.J)/self.length
        #     bend6Z = (6 * self.material.E * self.crossSection.Iz) / (self.length ** 2 * (1 + Constant_y))
        #     bend6Y = (6 * self.material.E * self.crossSection.Iy) / (self.length ** 2 * (1 + Constant_z))
        #     bendAlt4Z = ((4 + Constant_y) * self.material.E * self.crossSection.Iz) / (self.length * (1 + Constant_y))
        #     bendAlt4Y = ((4 + Constant_z)*self.material.E * self.crossSection.Iy) / (self.length * (1 + Constant_z))
        #     bendAlt2Z = ((2 - Constant_y) * self.material.E * self.crossSection.Iz) / (self.length * (1 + Constant_y))
        #     bendAlt2Y = ((2 - Constant_z) * self.material.E * self.crossSection.Iy) / (self.length * (1 + Constant_z))
        #
        #
        #     KE = np.zeros(12,12)
        #
        #     KE[0][0] = T
        #     KE[6][0] = -T
        #     KE[0][6] = -T
        #     KE[6][6] = T
        #
        #     KE[1][1] = bend12Z
        #     KE[2][2] = bend12Y
        #     KE[3][3] = twisting
        #     KE[1][5] = bend6Z
        #     KE[2][4] = - bend6Y
        #     KE[4][2] = - bend6Y
        #     KE[5][1] = bend6Z
        #     KE[4][4] = bendAlt4Y
        #     KE[5][5] = bendAlt4Z
        #
        #     KE[1][7] = - bend12Z
        #     KE[2][8] = - bend12Y
        #     KE[3][9] = - twisting
        #     KE[1][11] = bend6Z
        #     KE[2][10] = - bend6Y
        #     KE[4][8] = bend6Y
        #     KE[5][7] = - bend6Z
        #     KE[4][10] = bendAlt2Y
        #     KE[5][11] = bendAlt2Z
        #
        #     KE[7][1] = - bend12Z
        #     KE[8][2] = - bend12Y
        #     KE[9][3] = - twisting
        #     KE[7][5] = - bend6Z
        #     KE[8][4] = bend6Y
        #     KE[10][2] = - bend6Y
        #     KE[11][1] = bend6Z
        #     KE[10][4] = bendAlt2Y
        #     KE[11][5] = bendAlt2Z
        #
        #     KE[7][7] = bend12Z
        #     KE[8][8] = bend12Y
        #     KE[9][9] = twisting
        #     KE[7][11] = - bend6Z
        #     KE[8][10] = bend6Y
        #     KE[10][8] = bend6Y
        #     KE[11][7] = - bend6Z
        #     KE[10][10] = bendAlt4Y
        #     KE[11][11] = bendAlt4Z









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