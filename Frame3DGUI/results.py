
class NodeDeflections:
    def __init__(self, DX, DY, DZ, RX, RY, RZ):
        self.DX = DX
        self.DY = DY
        self.DZ = DZ
        self.RX = RX
        self.RY = RY
        self.RZ = RZ

class Reactions:
    def __init__(self, reactions):
        self.RX = reactions[:,0]
        self.RY = reactions[:,1]
        self.RZ = reactions[:,2]
        self.MX = reactions[:,3]
        self.MY = reactions[:,4]
        self.MZ = reactions[:,5]

class MaximumInternalForces:
    def __init__(self, internalForces):
        self.internalForces = internalForces
        F, M = internalForces
        FX, FY, FZ = F
        MX, MY, MZ = M

        self.FX = FX[0]
        self.FX_case = FX[1]
        self.FY = FY[0]
        self.FY_case = FY[1]
        self.FZ = FZ[0]
        self.FZ_case = FZ[1]

        self.MX = MX[0]
        self.MX_case = MX[1]
        self.MY = MY[0]
        self.MY_case = MY[1]
        self.MZ = MZ[0]
        self.MZ_case = MZ[1]

class Results:
    def __init__(self):
        self.nodalDeflections = []
        self.weight = None
        self.overallWeight = None
        self.reactions = []
        self.maxInternalForces = None

    def addNodalDeflections(self, DX, DY, DZ, RX, RY, RZ):
        for case_id in range(len(DX)):
            self.nodalDeflections.append(NodeDeflections(DX[case_id], DY[case_id], DZ[case_id],
                                                         RX[case_id], RY[case_id], RZ[case_id]))

    def addWeight(self, weight):
        self.weight = weight
        self.overallWeight = sum(weight)

    def addReactions(self, reactions):
        for case_id in range(len(reactions)):
            self.reactions.append(Reactions(reactions[case_id]))

    def addInternalForces(self, internalForces):
        self.maxInternalForces = MaximumInternalForces(internalForces)