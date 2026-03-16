
class Node:
    def __init__(self, cords, index):

        self.cords = cords
        self.index = index

        self.numDimensions = None
        self.calcNumDimensions()

        self.maxDisplacement = None
        self.displacememnt = None

        self.unitDisplacement = None

        self.activeDisplacementConstraint = None
        self.AmountViolated = None

    def calcNumDimensions(self):
        self.numDimensions = len(self.cords)

    def CheckActiveDisplacementConstraint(self):
        self.activeDisplacementConstraint = []
        self.AmountViolated = []

        for i in range(len(self.displacememnt)):
            self.activeDisplacementConstraint.append([])
            self.AmountViolated.append([])
            for j in range(len(self.cords)):
                if self.maxDisplacement[j] < abs(self.displacememnt[i][j]):
                    self.activeDisplacementConstraint[i].append(True)
                else:
                    self.activeDisplacementConstraint[i].append(False)

                if self.activeDisplacementConstraint[i][j]:
                    if abs(self.displacememnt[i][j]) == self.displacememnt[i][j]:
                        self.AmountViolated[i].append(self.displacememnt[i] - self.maxDisplacement)
                    else:
                        self.AmountViolated[i].append(-1 *  self.displacememnt[i] - self.maxDisplacement)