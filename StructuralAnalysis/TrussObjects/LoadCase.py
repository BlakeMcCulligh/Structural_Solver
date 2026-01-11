import numpy as np


class LoadCase:
    def __init__(self):

        self.nodePointLoads = []
        self.applyedForce = None

        self.internalForces = None
        self.deformation = None
        self.reactions = None

    def addNodePointLoad(self, nodeIndex, direction, magnitued):
        """

        :param nodeIndex: The index of the node the load is applyed to
        :param direction: (x = 1, y = 2, z = 3)
        :param magnitued: the magnitued of the load being applyed. A negative value indecated in the negative direcction of the direction indicated.
        :return:
        """

        self.nodePointLoads.append([nodeIndex, direction, magnitued])

    def arangeForces(self, numDimensions, numNodes):

        self.applyedForce = np.zeros(shape=(numDimensions * numNodes, 1))

        for i in range(len(self.nodePointLoads)):
            self.applyedForce[numDimensions * (self.nodePointLoads[i][0]) + self.nodePointLoads[i][1] - 1] = self.nodePointLoads[i][2]