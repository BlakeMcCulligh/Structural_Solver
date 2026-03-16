from math import pi

import numpy as np

from StructuralAnalysis.OLD.HelperFuctions import getQubicRoots

class Member:
    def __init__(self, nodes, crossSection, material):
        self.nodes = nodes
        self.crossSection = crossSection
        self.material = material
        self.length = None
        self.GlobalStiffnessMatrix = None

        self.calcLength()
        self.assembaleStifnessMatrix()

        self.normalForces = None
        self.normalStress = None

        self.comprtionResistance = None
        self.tentionResistance = None

        self.displacementVectors = []
        self.unitDisplacementVectors = []

    def calcLength(self):
        length = 0
        for i in range(self.nodes[0].numDimensions):
            length += (self.nodes[1].cords[i] - self.nodes[0].cords[i]) ** 2

        self.length = length ** 0.5

    def assembaleStifnessMatrix(self):

        if self.nodes[0].numDimensions == 3:
            # local Stiffness Matrix
            K_mPrime = np.zeros((6,6))

            K_mPrime[0][0] = 1
            K_mPrime[3][0] = -1
            K_mPrime[0][3] = -1
            K_mPrime[3][3] = 1

            mult = self.material.E * self.crossSection.A / self.length
            K_mPrime = mult * K_mPrime


            # Transformation matrix
            Cx = (self.nodes[1].cords[0] - self.nodes[0].cords[0]) / self.length
            Cy = (self.nodes[1].cords[1] - self.nodes[0].cords[1]) / self.length
            Cz = (self.nodes[1].cords[2] - self.nodes[0].cords[2]) / self.length

            T_m = np.zeros((6,6))

            T_m[0][0] = Cx
            T_m[3][3] = Cx

            T_m[0][1] = Cy
            T_m[3][4] = Cy

            T_m[0][2] = Cz
            T_m[3][5] = Cz

            # Convert to global stiffness matrix
            T_m_T = T_m.transpose()
            K_m = T_m_T @ K_mPrime @ T_m

            self.GlobalStiffnessMatrix = K_m

        else:
            # local Stiffness Matrix
            K_mPrime = np.zeros((4,4))

            K_mPrime[0][0] = 1
            K_mPrime[2][0] = -1
            K_mPrime[0][2] = -1
            K_mPrime[2][2] = 1

            mult = self.material.E * self.crossSection.A / self.length
            K_mPrime = mult * K_mPrime

            # Transformation matrix
            C = (self.nodes[1].cords[0] - self.nodes[0].cords[0]) / self.length
            S = (self.nodes[1].cords[1] - self.nodes[0].cords[1]) / self.length

            T_m = np.zeros((4,4))
            T_m[0][0] = C
            T_m[1][1] = C
            T_m[2][2] = C
            T_m[3][3] = C

            T_m[0][1] = -S
            T_m[1][0] = S
            T_m[2][3] = -S
            T_m[3][2] = S

            # Convert to global stiffness matrix
            T_m_T = T_m.transpose()
            K_m = T_m @ K_mPrime @ T_m_T

            self.GlobalStiffnessMatrix = K_m

    def calcNormalStress(self):
        self.normalStress = []
        for nForce in self.normalForces:
            self.normalStress.append(nForce / self.crossSection.A)

    def calcTensionResistance(self):

        self.tentionResistance = 0.9 * self.crossSection.Ag * self.material.Fy

    def calcComprestionResistance(self):
        Fe = self.getFe()
        Lambda = (self.material.Fy / Fe) ** 0.5
        self.comprtionResistance = 0.9 * self.crossSection.A * self.material.Fy / (1 + Lambda ** (2 * self.material.n)) ** (
                    1 / self.material.n)

    def getFe(self):
        ro2 = self.crossSection.xo ** 2 + self.crossSection.yo ** 2 + self.crossSection.rx ** 2 + self.crossSection.ry ** 2

        Fex = (pi ** 2) * self.material.E / (self.length / self.crossSection.rx)
        Fey = (pi ** 2) * self.material.E / (self.length / self.crossSection.ry)
        Fez = ((pi ** 2 * self.material.E * self.crossSection.Cw) / (
                    self.length) ** 2 + self.material.G * self.crossSection.J) * 1 / self.crossSection.A * ro2

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

    def calcDisplacememntVectors(self):
        self.displacementVectors = []
        for i in range(len(self.nodes[0].displacememnt)):
            d = self.nodes[0].displacement[i]
            d.extend(self.nodes[1].displacement[i])
            self.displacementVectors.append(d)

    def calcUnitDisplacementVector(self):
        self.UnitDisplacementVectors = []
        for i in range(len(self.nodes[0].unitDisplacement)):
            d = self.nodes[0].unitDisplacement[i]
            d.extend(self.nodes[1].unitDisplacement[i])
            self.UnitDisplacementVectors.append(d)





