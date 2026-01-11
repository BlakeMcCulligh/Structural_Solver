import numpy as np


class Zone:
    def __init__(self, nodes, edges):
        self.nodes = np.array(nodes)
        self.edges = np.array(edges)

    def getMostExtreamPoints(self):
        minX = self.nodes[:, 0].min()
        maxX = self.nodes[:, 0].max()
        minY = self.nodes[:, 1].min()
        maxY = self.nodes[:, 1].max()

        return minX, minY, maxX, maxY






