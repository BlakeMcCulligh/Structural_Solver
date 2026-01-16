import math

import numpy as np
from matplotlib import pyplot as plt

from TopologyOptimization.Zone import Zone

class Mesh:
    def __init__(self, zone, maxSpassing, CollinerarTollerance):
        self.zone = zone
        self.maxSpassing = maxSpassing
        self.nodes = None
        self.members = None

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

        # for i in range(len(self.nodes)):
        #     for j in range(len(self.nodes)):
        #         if i != j and self.nodes[i][0] == self.nodes[j][0] and self.nodes[i][1] == self.nodes[j][1]:
        #             print("node Deleted")
        #             del self.nodes[i]

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

    def generateLevel1Members(self):
        A = np.zeros((len(self.nodes),len(self.nodes)))

        for i in range(len(self.nodes)):
            x = self.nodes[i][0]
            y = self.nodes[i][1]
            closestXIndex = None
            disToClosestX = math.inf
            closestYIndex = None
            disToClosestY = math.inf
            for j in range(len(self.nodes)):

                if i != j:
                    if x == self.nodes[j][0]:
                        dist = y - self.nodes[j][1]
                        if dist < 0:
                            if abs(dist) < disToClosestX:
                                disToClosestX = dist
                                closestXIndex = j

                    elif y == self.nodes[j][1]:
                        dist = x - self.nodes[j][0]
                        if dist < 0:
                            if abs(dist) < disToClosestY:
                                disToClosestY = dist
                                closestYIndex = j

            if closestXIndex is not None:
                A[i,closestXIndex] = 1
            if closestYIndex is not None:
                A[i, closestYIndex] = 1

        A = A + A.T
        n = self.getListOfUnderseportedNodes(A)
        if len(n) != 0:
            repet = True
        else:
            repet = False
        while repet:

            node = self.nodes[n[0]]

            distClosest = math.inf
            closestIndex = None
            for i in range(len(self.nodes)):
                dist = ((node[0]-self.nodes[i][0])**2 + (node[1]-self.nodes[i][1])**2)**0.5
                if dist < distClosest and n[0] != i and A[n[0],i] == 0:
                    distClosest = dist
                    closestIndex = i

            A[n[0],closestIndex] = 1

            A = A + A.T
            n = self.getListOfUnderseportedNodes(A)
            if len(n) == 0:
                repet = False

        self.H = A
        self.A1 = A
        self.convertMatrixToMembers(A)

    def convertMatrixToMembers(self, A):
        self.members = []

        for i in range(len(A)):
            for j in range(len(A[i])):
                if A[i][j] > 0 and i > j:
                    self.members.append([i,j])

    def getListOfUnderseportedNodes(self, A):
        n = [] #TODO if corner node make only need 2 members
        for i in range(len(A)):
            numMembers = 0
            for j in range(len(A[i])):
                if A[i][j] > 0 and i != j:
                    numMembers += 1
            if numMembers < 3:
                n.append(i)

        return n

    def generateUpperLevelMembers(self, level):
        An = self.A1
        Anlow = None
        for i in range(level - 1):
            Anlow = An
            An = An @ self.A1

        Gn = An - Anlow

        #TODO remove all new members that are within a tolerance of colinarity of an already existing members
        # figure out what is happening at the top of the mesh.
        #self.members = []
        # m = np.array(self.members)
        # n = np.array(self.nodes)
        #
        # d = [n[m[:, 1],0] - n[m[:, 0],0], n[m[:, 1],1] - n[m[:, 0],1]]
        # l = ((d[0]) ** 2 + (d[1]) ** 2) ** 0.5
        #
        # unitVectors = np.array(d)/l
        #
        # for i in range(len(Gn)):
        #     for j in range(len(Gn)):
        #         if i != j and Gn[i,j] > 0:
        #             dnew = [n[i,0] - n[j,0], n[i,1] - n[j,1]]
        #             lnew = ((dnew[0]) ** 2 + (dnew[1]) ** 2) ** 0.5
        #             unitVectorNew = np.array(dnew) / lnew
        #
        #             unitVectorsNew = np.array([unitVectorNew] * len(unitVectors.T))
        #
        #             C = unitVectors * unitVectorsNew.T
        #
        #             #print(np.max(C))
        #             if np.max(C) > self.CollinerarTollerance:
        #                 Gn[i,j] = 0

        self.H = self.H + Gn
        self.H = self.H + self.H.T

        self.convertMatrixToMembers(self.H)




nodes = [[0,0],[0,100],[100,100],[100,0]]
edges = [[0,1],[1,2],[2,3],[3,0]]

Zone = Zone(nodes, edges)

MaxSpassing = [5,5]

mesh = Mesh(Zone, MaxSpassing, 0.999)
mesh.generateNodes()
mesh.generateLevel1Members()
# mesh.generateUpperLevelMembers(2)
# mesh.generateUpperLevelMembers(3)
# mesh.generateUpperLevelMembers(4)
mesh.print()