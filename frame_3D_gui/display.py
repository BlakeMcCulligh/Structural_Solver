"""
Handels the storag of things to be printed and the converting of said objects to a
printeble type: (node, line, and tris)
"""

import numpy as np

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

class Display:
    """
    Holds all things that are to be printed.
    """

    def __init__(self):
        """
        Display Object Constructor.
        """

        self.scale = [0.1]
        self.res = [8]

        self.Nodes = [] # shape: (# nodes, 3)
        self.Members = [] # shape: (# members, 6)

        self.supports = []  # shape: (# supports, 2: [[x,y,z],[dx,dy,dz,rx,ry,rz]])
        self.releces = []  # shape: (# releces, 2: [[x1,y1,z2,x2,y2,z2],12])
        self.point_loads = []  # shape: (# loads, 2: [[x,y,z],[dx,dy,dz,rx,ry,rz]])
        self.dist_loads = []  # shape: (# loads x directions, 2:[[x1,y1,z2,x2,y2,z2], [# seg, val]])

        self.node_deflections = []  # shape: (# delfelctions, 2: [[x,y,z],[dx,dy,dz,rx,ry,rz]])
        self.member_deflections = []  # shape: (# members, # seg, 3: [x,y,z]])
        self.reactions = []  # shape: (# reactions, 2: [[x,y,z],[dx,dy,dz,rx,ry,rz]])
        self.internal_loads = []  # shape: (# loads x directions, 2:[[x1,y1,z2,x2,y2,z2], [# seg, val]])

        self.text = [] # shape: (# text, 7: [string, x1,y1,z1, x2,y2,z2])

        self.PrintNodes: np.ndarray = np.empty((0, 3))
        self.PrintLines: np.ndarray = np.empty((0, 2, 3))
        self.PrintSurfaceTri: np.ndarray = np.empty((0, 3, 3))
        self.PrintSolidTri: np.ndarray = np.empty((0, 3, 3))

    def AddNode(self, node):
        """
        Addes a node to be displeyed in the 3D rendering.

        :param node: list. Node to be addes. shape: (3)
        """
        self.Nodes.append(node)
        self.ConvertToPrint()

    def AddMember(self, list_nodes, member):
        """
        Addes a member to be displeyed in the 3D rendering.

        :param list_nodes: list. List of all nodes in the frame. shape: (3, # nodes)
        :param member: list. Line end Node indeces in the print Node array [i Node, j Node].
        """

        nodes = [list_nodes[0][member[0]], list_nodes[1][member[0]], list_nodes[2][member[0]],
                 list_nodes[0][member[1]], list_nodes[1][member[1]], list_nodes[2][member[1]]]

        self.Members.append(nodes)
        self.ConvertToPrint()

    def AddSupports(self, list_nodes, support):
        location = [list_nodes[0][int(support[0])], list_nodes[1][int(support[0])], list_nodes[2][int(support[0])]]
        supports = support[1:]
        self.supports.append([location, supports])
        self.ConvertToPrint()

    def ConvertToPrint(self):
        """
        Converts all arrays of things to be displayed into the individual parts that can be printeded.
        """

        self.PrintNodes: np.ndarray = np.empty((0, 3))
        self.PrintLines: np.ndarray = np.empty((0, 2, 3))
        self.PrintSurfaceTri: np.ndarray = np.empty((0, 3, 3))
        self.PrintSolidTri: np.ndarray = np.empty((0, 3, 3))

        # nodes
        for i in range(len(self.Nodes)):
            self.PrintNodes = np.vstack((self.PrintNodes, self.Nodes[i]))

        # members
        for i in range(len(self.Members)):
            self.PrintLines = np.vstack((self.PrintLines,
                                         [[[self.Members[i][0], self.Members[i][1], self.Members[i][2]],
                                          [self.Members[i][3], self.Members[i][4], self.Members[i][5]]]]))

        # supports
        self._convert_supports()

    def _convert_supports(self):
        s = self.scale[0]



        X_circle = [s / 4 * np.cos(2 * np.pi * i / self.res[0]) for i in range(self.res[0])]
        Y_circle = [s / 4 * np.sin(2 * np.pi * i / self.res[0]) for i in range(self.res[0])]

        for support in self.supports:
            x = support[0][0]
            y = support[0][1]
            z = support[0][2]

            if all(support[1]): # fully fixed
                # solid horizontal rectangel at top
                n1 = [x - s, y - s, z]
                n2 = [x - s, y + s, z]
                n3 = [x + s, y + s, z]
                n4 = [x + s, y - s, z]
                tri = [[n1, n2, n3], [n1, n4, n3]]
                for i in range(len(tri)): self.PrintSurfaceTri = np.vstack((self.PrintSurfaceTri, [tri[i]]))


            else:
                if support[1][2]:  # virtical support

                    if support[1][0]:  # X horizontal support
                        if support[1][4]:  # Y rotaition support
                            # x axis rectangle
                            n1 = [x - s, y, z]
                            n2 = [x - s, y, z - s]
                            n3 = [x + s, y, z - s]
                            n4 = [x + s, y, z]
                            self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n3, n4]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n4, n1]]))

                        else:
                            # x axis triangle
                            n1 = [x, y, z]
                            n2 = [x - s, y, z - s]
                            n3 = [x + s, y, z - s]
                            self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n3, n1]]))

                    else:
                        # circle and under-line
                        X = [x + X_circle[i] for i in range(self.res[0])]
                        Y = [y] * self.res[0]
                        Z = [z - s*0.75 + Y_circle[i] for i in range(self.res[0])]
                        for i in range(len(X) - 1):
                            self.PrintLines = np.vstack((self.PrintLines, [[[X[i], Y[i], Z[i]],
                                                                           [X[i + 1], Y[i + 1], Z[i + 1]]]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[[X[len(X) - 1], Y[len(X) - 1], Z[len(X) - 1]],
                                                                        [X[0], Y[0], Z[0]]]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[[x-s/4,y,z-s],[x+s/4,y,z-s]]]))

                        if support[1][4]:  # Y rotaition support
                            # x-axis half rectangle
                            n1 = [x - s/4, y, z]
                            n2 = [x - s/4, y, z - s/2]
                            n3 = [x + s/4, y, z - s/2]
                            n4 = [x + s/4, y, z]
                            self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n3, n4]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n4, n1]]))

                        else:
                            # x-axis half triangle
                            n1 = [x, y, z]
                            n2 = [x - s/4, y, z - s/2]
                            n3 = [x + s/4, y, z - s/2]
                            self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n3, n1]]))

                    if support[1][1]:  # y horizontal support
                        if support[1][3]:  # X rotaition support
                            # y axis rectangle
                            n1 = [x, y - s, z]
                            n2 = [x, y - s, z - s]
                            n3 = [x, y + s, z - s]
                            n4 = [x, y + s, z]
                            self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n3, n4]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n4, n1]]))
                        else:
                            # y axis triangle
                            n1 = [x, y, z]
                            n2 = [x, y - s, z - s]
                            n3 = [x, y + s, z - s]
                            self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n3, n1]]))

                    else:
                        # circle and under-line
                        X = [x] * self.res[0]
                        Y = [y + X_circle[i] for i in range(self.res[0])]
                        Z = [z - s*0.75 + Y_circle[i] for i in range(self.res[0])]
                        for i in range(len(X) - 1):
                            self.PrintLines = np.vstack((self.PrintLines, [[[X[i], Y[i], Z[i]],
                                                                           [X[i + 1], Y[i + 1], Z[i + 1]]]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[[X[len(X) - 1], Y[len(X) - 1], Z[len(X) - 1]],
                                                                        [X[0], Y[0], Z[0]]]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[[x, y - s/4, z - s], [x, y + s/4, z - s]]]))

                        if support[1][3]:  # X rotaition support
                            # y-axis half rectangle
                            n1 = [x, y - s/4, z]
                            n2 = [x, y - s/4, z - s / 2]
                            n3 = [x, y + s/4, z - s / 2]
                            n4 = [x, y + s/4, z]
                            self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n3, n4]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n4, n1]]))
                        else:
                            # y-axis half triangle
                            n1 = [x, y, z]
                            n2 = [x, y - s/4, z - s / 2]
                            n3 = [x, y + s/4, z - s / 2]
                            self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n3, n1]]))

                    if support[1][5]:  # Z rotaition support
                        # horizontal rectangel under
                        n1 = [x - s, y - s, z - s]
                        n2 = [x - s, y + s, z - s]
                        n3 = [x + s, y + s, z - s]
                        n4 = [x + s, y - s, z - s]
                        self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n3, n4]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n4, n1]]))

                else:

                    if support[1][0]:  # x horizontal support
                        # circle and under-line
                        X = [x - s * 0.75 + X_circle[i] for i in range(self.res[0])]
                        Y = [y + Y_circle[i] for i in range(self.res[0])]
                        Z = [z] * self.res[0]
                        for i in range(len(X) - 1):
                            self.PrintLines = np.vstack((self.PrintLines, [[[X[i], Y[i], Z[i]],
                                                                           [X[i + 1], Y[i + 1], Z[i + 1]]]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[[X[len(X) - 1], Y[len(X) - 1], Z[len(X) - 1]],
                                                                        [X[0], Y[0], Z[0]]]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[[x - s, y - s/4, z], [x - s, y + s/4, z]]]))

                        if support[1][5]:  # Z rotaition support
                            # horizontal half rectangel
                            n1 = [x - s/2, y - s/4, z]
                            n2 = [x - s/2, y + s/4, z]
                            n3 = [x, y + s/4, z]
                            n4 = [x, y - s/4, z]
                            self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n3, n4]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n4, n1]]))

                        else:
                            # horizontal half triangle
                            n1 = [x - s / 2, y - s/4, z]
                            n2 = [x - s / 2, y + s/4, z]
                            n3 = [x, y, z]
                            self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n3, n1]]))


                    if support[1][1]:
                        # circle and under-line
                        X = [x + X_circle[i] for i in range(self.res[0])]
                        Y = [y - s * 0.75 + Y_circle[i] for i in range(self.res[0])]
                        Z = [z] * self.res[0]
                        for i in range(len(X) - 1):
                            self.PrintLines = np.vstack((self.PrintLines, [[[X[i], Y[i], Z[i]],
                                                                           [X[i + 1], Y[i + 1], Z[i + 1]]]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[[X[len(X) - 1], Y[len(X) - 1], Z[len(X) - 1]],
                                                                        [X[0], Y[0], Z[0]]]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[[x - s/4, y - s, z], [x + s/4, y - s, z]]]))

                        if support[1][5]:  # Z rotaition support
                            # horizontal half rectangel
                            n1 = [x - s/4, y - s / 2, z]
                            n2 = [x + s/4, y - s / 2, z]
                            n3 = [x + s/4, y, z]
                            n4 = [x - s/4, y, z]
                            self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n3, n4]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n4, n1]]))
                        else:
                            # horizontal half triangle
                            n1 = [x - s/4, y - s / 2, z]
                            n2 = [x + s/4, y - s / 2, z]
                            n3 = [x, y, z]
                            self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n3, n1]]))

                    if support[1][3]:  # X rotaition support
                        # y-axis half rectangle
                        n1 = [x, y - s, z]
                        n2 = [x, y - s, z - s / 2]
                        n3 = [x, y + s, z - s / 2]
                        n4 = [x, y + s, z]
                        self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n3, n4]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n4, n1]]))

                        # y-axis half rectangle
                        n1 = [x, y - s, z]
                        n2 = [x, y - s, z - s / 2]
                        n3 = [x, y + s, z - s / 2]
                        n4 = [x, y + s, z]
                        self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n3, n4]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n4, n1]]))

                    if support[1][4]:  # Y rotaition support
                        # circle and under-line
                        X = [x + X_circle[i] for i in range(self.res[0])]
                        Y = [y] * self.res[0]
                        Z = [z - s * 0.75 + Y_circle[i] for i in range(self.res[0])]
                        for i in range(len(X) - 1):
                            self.PrintLines = np.vstack((self.PrintLines, [[[X[i], Y[i], Z[i]],
                                                                           [X[i + 1], Y[i + 1], Z[i + 1]]]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[[X[len(X) - 1], Y[len(X) - 1], Z[len(X) - 1]],
                                                                        [X[0], Y[0], Z[0]]]]))

                        # x-axis half rectangle
                        n1 = [x - s/4, y, z]
                        n2 = [x - s/4, y, z - s / 2]
                        n3 = [x + s/4, y, z - s / 2]
                        n4 = [x + s/4, y, z]
                        self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n3, n4]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n4, n1]]))

                    if support[1][5] and not (support[1][0] and support[1][1]):  # Z rotaition support
                        # horizontal rectangel
                        n1 = [x - s, y - s, z]
                        n2 = [x - s, y + s, z]
                        n3 = [x + s, y + s, z]
                        n4 = [x + s, y - s, z]
                        self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n3, n4]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n4, n1]]))

    # todo add inputs for arrays
    # todo add converter to convert everything to the print arrays
    # todo rework all the adders
    # todo popup to select what should be printed