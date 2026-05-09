"""
Holds the object the stores all the results, incuding sub objects.
"""

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

class NodeDeflections:
    """
    Node deflections object.
    """

    def __init__(self, DX, DY, DZ, RX, RY, RZ):
        """
        Node deflection object constructer.

        :param DX: X direction deflections. Shape: (# nodes)
        :param DY: Y direction deflections. Shape: (# nodes)
        :param DZ: Z direction deflections. Shape: (# nodes)
        :param RX: X direction rotation deflections. Shape: (# nodes)
        :param RY: Y direction rotation deflections. Shape: (# nodes)
        :param RZ: Z direction rotation deflections. Shape: (# nodes)
        """

        self.DX = DX
        self.DY = DY
        self.DZ = DZ
        self.RX = RX
        self.RY = RY
        self.RZ = RZ

class Reactions:
    """
    Reactions object.
    """

    def __init__(self, reactions):
        """
        Reactions object constructer.

        :param reactions: Node ractions. Shape: (# Nodes, 6)
        """

        self.RX = reactions[:,0]
        self.RY = reactions[:,1]
        self.RZ = reactions[:,2]
        self.MX = reactions[:,3]
        self.MY = reactions[:,4]
        self.MZ = reactions[:,5]

class MaximumInternalForces:
    """
    Maximum internal forces for each member object.
    """

    def __init__(self, internalForces):
        """
        Maximum internal forces for each member object constructor.

        :param internalForces: Maximum forces in each member. Shape: (2,3,2)
        """

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
    """
    Frame structural analysis results object.
    """

    def __init__(self):
        """
        Results object constructor.
        """

        self.NodalDeflections = []
        self.Weight = None
        self.OverallWeight = None
        self.Reactions = []
        self.MaxInternalForces = None

    def AddNodalDeflections(self, DX, DY, DZ, RX, RY, RZ):
        """
        Adds node deflections to the results object.

        :param DX: X direction deflections. Shape: (# casses, # nodes)
        :param DY: Y direction deflections. Shape: (# casses, # nodes)
        :param DZ: Z direction deflections. Shape: (# casses, # nodes)
        :param RX: X direction rotation deflections. Shape: (# casses, # nodes)
        :param RY: Y direction rotation deflections. Shape: (# casses, # nodes)
        :param RZ: Z direction rotation deflections. Shape: (# casses, # nodes)
        """

        for case_id in range(len(DX)):
            self.NodalDeflections.append(NodeDeflections(DX[case_id], DY[case_id], DZ[case_id],
                                                         RX[case_id], RY[case_id], RZ[case_id]))

    def AddWeight(self, weight):
        """
        Adds weigits to the results object.

        :param weight: List of wights of each member.
        """

        self.Weight = weight
        self.OverallWeight = sum(weight)

    def AddReactions(self, reactions):
        """
        Adds reactions to the results object.

        :param reactions: Node ractions. Shape: (# cases, # Nodes, 6)
        """

        for case_id in range(len(reactions)):
            self.Reactions.append(Reactions(reactions[case_id]))

    def AddInternalForces(self, internalForces):
        """
        Adds the maximum forces in each momber to the results object.

        :param internalForces: Maximum forces in each member. Shape: (2,3,2)
        """

        self.MaxInternalForces = MaximumInternalForces(internalForces)