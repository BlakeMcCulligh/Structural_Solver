import math

import numpy as np
from matplotlib import pyplot as plt

from TopologyOptimization.Optimize import TrussTopologyOptimization
from TopologyOptimization.Zone import Zone

class Mesh:
    def __init__(self, zone, maxSpassing, CollinerarTollerance):
        self.zone = zone
        self.maxSpassing = maxSpassing
        self.nodes = None
        self.members = []
        self.memberLengths = []

        self.CollinerarTollerance = CollinerarTollerance

        self.H = None
        self.A1 = None

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

        nodesToRemove = list(nodesToRemove)

        indecesToRemove = []
        for i in range(len(nodesToRemove)):
            for j in range(len(self.zone.nodes)):
                n1 = self.nodes[nodesToRemove[i]]
                n2 = self.zone.nodes[j]
                if n1[0] == n2[0] and n1[1] == n2[1]:
                    indecesToRemove.append(i)

        for index in sorted(indecesToRemove, reverse=True):
            del nodesToRemove[index]

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

    def print(self):
        plt.plot(self.nodes[:,0], self.nodes[:,1], 'o')

        for i in range(len(self.members)):
             xi, xf, = self.nodes[self.members[i][0], 0], self.nodes[self.members[i][1], 0]
             yi, yf, = self.nodes[self.members[i][0], 1], self.nodes[self.members[i][1], 1]
             plt.plot([xi, xf], [yi, yf], color='gray')

        plt.show()

    def addMembers(self, maxDistance):

        for i in range(len(self.nodes)):
            node = self.nodes[i]

            distance = ((self.nodes[:, 0] - node[0]) ** 2 + (self.nodes[:, 1] - node[1]) ** 2) ** 0.5
            index = []
            for j in range(len(self.nodes)):
                index.append(j)

            distance = distance.tolist()

            indexSorted = [x for _, x in sorted(zip(distance, index ))]

            nodesToAdd = []

            for ind in indexSorted:
                if i != ind and distance[ind] <= maxDistance:

                    sk = False
                    for j in range(len(self.members)):
                        if i == self.members[j][0] and ind == self.members[j][1]:
                            sk = True
                        elif i == self.members[j][1] and ind == self.members[j][0]:
                            sk = True

                    if not sk:
                        dontAdd = False
                        for j in range(len(nodesToAdd)):
                            if is_collinear(node, self.nodes[ind], self.nodes[nodesToAdd[j]]):
                                dontAdd = True
                        if not dontAdd:
                            nodesToAdd.append(ind)

            for node2Index in nodesToAdd:
                #TODO check if member goes outside of the bounds
                self.members.append([i, node2Index])

    def CalcMemberLengths(self):
        m = np.array(self.members)

        self.memberLengths = ((self.nodes[m[:, 1], 0] - self.nodes[m[:, 0], 0]) ** 2 + (self.nodes[m[:, 1], 1] - self.nodes[m[:, 0], 1]) ** 2 ) ** 0.5

def is_collinear(a, b, c):

    a = np.append(a, 0)
    b = np.append(b, 0)
    c = np.append(c, 0)

    AB = b - a
    AC = c - a

    cp = np.cross(AB, AC)

    return np.allclose(cp, 0)

nodes = [[0,0],[0,50],[50,50],[50,0]]
edges = [[0,1],[1,2],[2,3],[3,0]]

Zone = Zone(nodes, edges)

MaxSpassing = [10,10]

mesh = Mesh(Zone, MaxSpassing, 0.999)
mesh.generateNodes()
mesh.addMembers(30)
mesh.CalcMemberLengths()

#mesh.print()

volume_bound = 500
E = 200000
opt = TrussTopologyOptimization(volume_bound, np.array(mesh.members), np.array(mesh.nodes), mesh.memberLengths, E)
opt.optimize()
opt.printResults()