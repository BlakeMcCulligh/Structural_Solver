import numpy as np


class OptimizeTrussCrossSections:
    def __init__(self, truss):
        self.T = truss

        self.A = np.array([])

        self.deformationStepSize = 0.5 # TODO updater needed for end of each iteration

        self.ITER = 0

    def optimize(self):
        self.ITER = 0
        self.mainLoop()

    def mainLoop(self):

        self.loop2()
        self.loop3()
        self.loop4()
        self.loop5()
        self.loop6()
        self.loop7()

        self.ITER += 1

    def loop2(self):
        self.T.AssembleGlobalStiffnessMatrix()
        self.T.solveStiffnessMatix()
        self.T.calcStresses()

        # TODO find the maximum scaling factor obtained from the initial values for the cross-sections
        # Percent of the beginging Area and the range of areas alowed

    def loop3(self):
        a =1
        # TODO the design variables are estimated by satisfying the displacement constraints only, assuming continuous variables, and a new scaling factor is computed

    def loop4(self):
        a =1
        # map to cross-sections

    def loop5(self):
        a =1
        # all the members whose compressive
        # stress constraints are violated are scaled by SFS's and the
        # remaining members are scaled by y

        #all the members are mapped again to the discrete standard cross-sections from
        # the AISC W sections database

        #repeat loop5 until all compression members satisfy the stress and buckling constraints

    def loop6(self):
        a =1
        # calculate the weight of the structure and
        # resize the members according to Equation 13

    def loop7(self):
        a =1

    def deformationScalling(self, p):
        """

        :param p: Iteration Number starting at 1
        :return:
        """
        #self.T.checkActiveDisplacementConstraints()

        multPart1 = [] # k , j
        for k in range(len(self.T.loadCases)):
            multPart1.append(0)
            for n in range(len(self.T.nodes)):
                if self.T.nodes[n].activeDisplacementConstraint:
                    multPart1[k] = multPart1[k] + self.T.weight / np.array(self.T.nodes[n].displacementVectors[k])
        multPart1 = np.array(multPart1)

        multPart2outer = []  # i x j x k
        for i in range(len(self.T.memberGroupings)):
            numerator = [] # j, K

            for j in range(len(self.T.nodes[0].numDimensions)):
                numerator.append([])
                for k in range(len(self.T.loadCases)):  #
                    numerator[j].append(np.zeros((self.T.loadCases, 1)))
                    for m in self.T.memberGroupings[i]:
                        numerator[j][k] = numerator[j][k] + ((np.array(self.T.members[m].unitDisplacementVectors[j])).T @
                                                       self.T.members[m].GlobalStiffnessMatrix @
                                                       np.array(self.T.members[m].displacementVectors[k]))

            L = 0
            for m in self.T.memberGroupings[i]:
                L += self.T.members[m].length

            denominator = p * self.A[i] * L

            multPart2 = [] # j x k
            for j in range(len(self.T.nodes[0].numDimensions)):
                multPart2.append([])
                for k in range(len(self.T.loadCases)):
                    multPart2[j].append(numerator[j][k]/denominator)
            multPart2outer.append(np.array(multPart2))

        multiplyerBase = []
        for i in range(len(self.T.memberGroupings)):
            multiplyerBase.append(max(multPart1 @ multPart2outer[i]))

        multiplyer = np.array(multiplyerBase) ** self.deformationStepSize

        self.A = self.A * multiplyer

















