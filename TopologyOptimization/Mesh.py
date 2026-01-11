import math

import numpy as np

from TopologyOptimization.Zone import Zone

class Mesh:
    def __init__(self, zone, maxSpassing):
        self.zone = zone
        self.maxSpassing = maxSpassing
        self.nodes = None

    def generateNodes(self):
        minX, minY, maxX, maxY = self.zone.getMostExtreamPoints()

        xLength = maxX - minX
        yLength = maxY - minY

        numGapsX = math.ceil(xLength/self.maxSpassing[0])
        numGapsY = math.ceil(yLength/self.maxSpassing[1])

        dx = xLength/numGapsX
        dy = yLength/numGapsY

        self.nodes = []
        for i in range(numGapsY + 1):
            for j in range(numGapsX + 1):
                self.nodes.append([minX + dx * j, minY + dy * i])

        self.removeNodesOutsideZone()

        self.addNodesToEadgeOfZoneWhereNeeded(minX, minY, numGapsX, numGapsY, dx, dy)

    def removeNodesOutsideZone(self):
        nodeArray = np.array(self.nodes)

        x = nodeArray[:, 0]
        y = nodeArray[:, 1]

        xp = self.zone.nodes[:, 0]
        yp = self.zone.nodes[:, 1]

        n = len(self.zone.nodes)
        inside = np.zeros(len(nodeArray), dtype=bool)

        # Loop over polygon edges
        for i in range(n):
            j = (i + 1) % n

            xi, yi = xp[i], yp[i]
            xj, yj = xp[j], yp[j]

            crosses = (yi > y) != (yj > y)

            if (yj - yi) != 0:
                x_intersect = xi + (y - yi) * (xj - xi) / (yj - yi)
            else:
                x_intersect = x - 1

            intersect = crosses & (x < x_intersect)

            inside ^= intersect

        nodesToRemove = np.where(~inside)[0]

        self.nodes = np.delete(nodeArray, nodesToRemove, axis=0)

    def addNodesToEadgeOfZoneWhereNeeded(self, minX, minY, numGapsX, numGapsY, dx, dy):

        lineSagmentEnds = [self.zone.nodes[self.zone.edges[:, 0], :], self.zone.nodes[self.zone.edges[:, 1], :]]

        vectors = [lineSagmentEnds[1][:, 0] - lineSagmentEnds[0][:, 0], lineSagmentEnds[1][:, 1] - lineSagmentEnds[0][:, 1]]

        Intervals = [[],[]]
        for i in range(numGapsY):
            Intervals[1].append(minY + dy * i)

        for i in range(numGapsX):
            Intervals[0].append(minX + dx * i)

        Intervals = [np.array(Intervals[0]), np.array(Intervals[1])]

        vectors = np.array(vectors)
        for i in range(len(Intervals[0])):
            for j in range(len(lineSagmentEnds[0])):

                if vectors[0][j] != 0:
                    K = (Intervals[0][i] - lineSagmentEnds[0][j, 0]) / vectors[0][j]
                else:
                    K = 2

                if 0 <= K < 1:
                    point = lineSagmentEnds[0][j,:] + vectors[:,j] * K

                    if not self.doesPointExist(point):

                        point = np.array([[point[0], point[1]]])
                        self.nodes = np.append(self.nodes, point, axis =0)

        for i in range(len(Intervals[1])):
            for j in range(len(lineSagmentEnds[0])):

                if vectors[1][j] != 0:
                    K = (Intervals[1][i] - lineSagmentEnds[0][j, 1]) / vectors[1][j]
                else:
                    K = 2

                if 0 <= K < 1:
                    point = lineSagmentEnds[0][j,:] + vectors[:,j] * K

                    if not self.doesPointExist(point):

                        point = np.array([[point[0], point[1]]])
                        self.nodes = np.append(self.nodes, point, axis =0)

    def doesPointExist(self, point):

        for node in self.nodes:
            if node[0] == point[0] and node[1] == point[1]:
                return True

        return False



nodes = [[0,0],[0,100],[100,100],[100,0]]
edges = [[0,1],[1,2],[2,3],[3,0]]

Zone = Zone(nodes, edges)

MaxSpassing = [5,5]

mesh = Mesh(Zone, MaxSpassing)
mesh.generateNodes()