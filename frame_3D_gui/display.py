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
        self.releces = []  # shape: (# releces, 2: [[x1,y1,z1,x2,y2,z2],12])
        self.point_loads = []  # shape: (# loads, 2: [[x,y,z],[dx,dy,dz,rx,ry,rz]])
        self.dist_loads = []  # shape: (# loads x directions, 2:[[x1,y1,z2,x2,y2,z2], [wx1, wx2, wy1, wy2, wz1, wz2]])

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
        """
        Adds supports to be displeyed in the 3D rendering.

        :param list_nodes: list. List of all nodes in the frame. shape: (3, # nodes)
        :param support: list. Support nodes and what DOF are supported. shape: [i Node, DX, DY, DZ, RX, RY, RZ]
        """

        location = [list_nodes[0][int(support[0])], list_nodes[1][int(support[0])], list_nodes[2][int(support[0])]]
        supports = support[1:]
        self.supports.append([location, supports])
        self.ConvertToPrint()

    def AddReleces(self, list_members, list_nodes, releces):
        """
        Adds releces to be displeyed in the 3D rendering.

        :param list_members: list. List of all members in the frame. shape: (2: [i Node, j Node], # nodes)
        :param list_nodes: list. List of all nodes in the frame. shape: (3, # nodes)
        :param releces: list, Releces for a member.
                        shape: [i member, iDX, iDY, iDZ, iRX, iRY, iRZ, jDX, jDY, jDZ, jRX, jRY, jRZ]
        """

        n1 = [list_nodes[0][list_members[0][releces[0]]],
              list_nodes[1][list_members[0][releces[0]]],
              list_nodes[2][list_members[0][releces[0]]]]

        n2 = [list_nodes[0][list_members[1][releces[0]]],
              list_nodes[1][list_members[1][releces[0]]],
              list_nodes[2][list_members[1][releces[0]]]]

        locations = n1 + n2
        direction = releces[1:]
        self.releces.append([locations, direction])
        self.ConvertToPrint()

    def AddNodeLoads(self, list_nodes, node_loads):
        """
        Adds node loads to be displeyed in the 3D rendering.

        :param list_nodes: list. List of all nodes in the frame. shape: (3, # nodes)
        :param node_loads: list. Node load. shape: [i Node, PX, PY, PZ, MX, MY, MZ]
        """

        locations = [list_nodes[0][node_loads[0]],list_nodes[1][node_loads[0]],list_nodes[2][node_loads[0]]]
        load = node_loads[1:]
        self.point_loads.append([locations, load])

    def AddMemberPointLoads(self, list_members, list_nodes, member_point_loads):
        """
        Adds member point loads to be displeyed in the 3D rendering.

        :param list_members: list. List of all members in the frame. shape: (2: [i Node, j Node], # nodes)
        :param list_nodes: list. List of all nodes in the frame. shape: (3, # nodes)
        :param member_point_loads: list. Local member point load. shape: [i Member, x, PX, PY, PZ, MX, MY, MZ]
        """

        location = get_loc_on_member(int(member_point_loads[0]), int(member_point_loads[1]), list_members, list_nodes)
        Trans = get_member_t(int(member_point_loads[0]), list_members, list_nodes)[:3, :3]
        P_global = Trans.T @ np.array([member_point_loads[2],member_point_loads[3],member_point_loads[4]]) @ Trans
        M_global = Trans.T @ np.array([member_point_loads[5],member_point_loads[6],member_point_loads[7]]) @ Trans
        self.point_loads.append([location, np.concatenate((P_global,M_global))])

    def AddMemberDistLoads(self, list_members, list_nodes, member_dist_loads):
        """
        Adds member distributed loads to be displeyed in the 3D rendering.

        :param list_members: list. List of all members in the frame. shape: (2: [i Node, j Node], # nodes)
        :param list_nodes: list. List of all nodes in the frame. shape: (3, # nodes)
        :param member_dist_loads: list. Local member distributed load.
                                  shape: [i Member, x1, x2, wx1, wx2, wy1, wy2, wz1, wz2]
        """

        location1 = get_loc_on_member(member_dist_loads[0], member_dist_loads[1], list_members, list_nodes)
        location2 = get_loc_on_member(member_dist_loads[0], member_dist_loads[2], list_members, list_nodes)
        location = location1 + location2
        Trans = get_member_t(member_dist_loads[0], list_members, list_nodes)[:3, :3]
        P1_global = Trans.T @ np.array([member_dist_loads[3], member_dist_loads[5], member_dist_loads[7]]) @ Trans
        P2_global = Trans.T @ np.array([member_dist_loads[4], member_dist_loads[6], member_dist_loads[8]]) @ Trans
        load = [P1_global[0],P2_global[0],P1_global[1],P2_global[1],P1_global[2],P2_global[2]]
        self.dist_loads.append([location, load])

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

    # todo add converter to convert everything to the print arrays
    # todo rework all the adders
    # todo popup to select what should be printed

def get_loc_on_member(i_member: int, x, list_members, list_nodes):
    """
    Finds the global coords of where x lands on a member.

    :param i_member: int, index of member.
    :param x: float. location along member.
    :param list_members: list. List of all members in the frame. shape: (2: [i Node, j Node], # nodes)
    :param list_nodes: list. List of all nodes in the frame. shape: (3, # nodes)
    :return: Location on member.
    """

    n1 = np.array([list_nodes[0][int(list_members[0][i_member])],
                   list_nodes[1][int(list_members[0][i_member])],
                   list_nodes[2][int(list_members[0][i_member])]])

    n2 = np.array([list_nodes[0][int(list_members[1][i_member])],
                   list_nodes[1][int(list_members[1][i_member])],
                   list_nodes[2][int(list_members[1][i_member])]])

    vector = n2 - n1
    unit_vector = vector / np.linalg.norm(vector)

    return  n1 + unit_vector * x

def get_member_t(i_member, list_members, list_nodes):
    """
    Builds an array of the transformation matrices for the member.

    :param i_member: int, index of member.
    :param list_members: list. List of all members in the frame. shape: (2: [i Node, j Node], # nodes)
    :param list_nodes: list. List of all nodes in the frame. shape: (3, # nodes)
    :return: transformation matrix.
    """

    # Get the global coordinates for the two ends
    n1 = np.array([list_nodes[0][int(list_members[0][i_member])],
                   list_nodes[1][int(list_members[0][i_member])],
                   list_nodes[2][int(list_members[0][i_member])]])

    n2 = np.array([list_nodes[0][int(list_members[1][i_member])],
                   list_nodes[1][int(list_members[1][i_member])],
                   list_nodes[2][int(list_members[1][i_member])]])

    i = list_members[0][i_member]
    j = list_members[1][i_member]
    Xi, Yi, Zi = n1
    Xj, Yj, Zj = n2

    # Calculate the length of the member
    L = np.linalg.norm(n2-n1)

    # Calculate the direction cosines for the local x-axis
    x = [(Xj - Xi) / L, (Yj - Yi) / L, (Zj - Zi) / L]

    # Vertical members
    if np.isclose(Xi, Xj) and np.isclose(Zi, Zj):

        if Yj > Yi:
            y = [-1, 0, 0]
            z = [0, 0, 1]
        else:
            y = [1, 0, 0]
            z = [0, 0, 1]

    # Horizontal members
    elif np.isclose(Yi, Yj):

        y = [0, 1, 0]
        z = np.cross(x, y)

        z = np.divide(z, (z[0] ** 2 + z[1] ** 2 + z[2] ** 2) ** 0.5)

    # Members neither vertical nor horizontal
    else:

        proj = [Xj - Xi, 0, Zj - Zi]

        if Yj > Yi:
            z = np.cross(proj, x)
        else:
            z = np.cross(x, proj)

        z = np.divide(z, (z[0] ** 2 + z[1] ** 2 + z[2] ** 2) ** 0.5)

        y = np.cross(z, x)
        y = np.divide(y, (y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5)

    # Create the direction cosines matrix
    dirCos = np.array([x, y, z])

    # Build the transformation matrix
    transMatrix = np.zeros((12, 12))
    transMatrix[0:3, 0:3] = dirCos
    transMatrix[3:6, 3:6] = dirCos
    transMatrix[6:9, 6:9] = dirCos
    transMatrix[9:12, 9:12] = dirCos

    return transMatrix