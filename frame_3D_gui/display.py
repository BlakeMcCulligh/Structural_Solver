"""
Handles the storage of things to be printed and the converting of said objects to a
printable type: (node, line, and tris)
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, List, Union
import numpy as np

if TYPE_CHECKING:
    from window import MainWindow

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

    def __init__(self, window: MainWindow):
        """
        Display Object Constructor.
        """

        self.window = window

        self.scale = [0.1,0.1,0.1,1]
        self.res = [8]

        self.Nodes: List[List[float]] = [] # shape: (# nodes, 3)
        self.Members: List[List[float]] = [] # shape: (# members, 6)

        # shape: (# supports, 2: [[x,y,z],[dx,dy,dz,rx,ry,rz]])
        self.supports: List[List[List[Union[float, bool]]]] = []
        # shape: (# releces, 2: [[x1,y1,z1,x2,y2,z2],12,member i])
        self.releces: List[List[Union[int,List[Union[float, bool]]]]] = []
        # shape: (# loads, 2: [[x,y,z],[dx,dy,dz,rx,ry,rz]])
        self.point_loads: List[List[List[float]]] = []
        # shape: (# loads x directions, 2:[[x1,y1,z2,x2,y2,z2], [wx1, wx2, wy1, wy2, wz1, wz2]])
        self.dist_loads: List[List[List[float]]] = []

        self.node_deflections: list = []  # shape: (# deflections, 2: [[x,y,z],[dx,dy,dz,rx,ry,rz]])
        self.member_deflections: list = []  # shape: (# members, # seg, 3: [x,y,z]])
        self.reactions: list = []  # shape: (# reactions, 2: [[x,y,z],[dx,dy,dz,rx,ry,rz]])
        self.internal_loads: list = []  # shape: (# loads x directions, 2:[[x1,y1,z2,x2,y2,z2], [# seg, val]])

        self.text: list = [] # shape: (# text, 7: [string, x1,y1,z1, x2,y2,z2])

        self.PrintNodes: np.ndarray = np.empty((0, 3))
        self.PrintLines: np.ndarray = np.empty((0, 2, 3))
        self.PrintSurfaceTri: np.ndarray = np.empty((0, 3, 3))
        self.PrintSolidTri: np.ndarray = np.empty((0, 3, 3))
        self.PrintText: list = [[],np.empty((0,6))] # shape: [[strings], (# text, 6)]

    def AddNode(self, node: List[float], i = None) -> None:
        """
        Adds a node to be displayed in the 3D rendering.

        :param node: list. Node to be added. shape: (3)
        :param i: index of the node. optional, only needed if to be displayed.
        """
        self.Nodes.append(node)

        if i is not None:
            n = np.array(node)
            n[2] += self.scale[3]
            cross = np.cross(self.window.LookDir, self.window.Up)

            n1 = (n + cross).tolist()
            n2 = (n - cross).tolist()
            self.text.append([i]+n1+n2)

        self.ConvertToPrint()

    def AddMember(self, list_nodes: List[List[float]], member: List[int]) -> None:
        """
        Adds a member to be displayed in the 3D rendering.

        :param list_nodes: list. List of all nodes in the frame. shape: (3, # nodes)
        :param member: list. Line end Node indeces in the print Node array [i Node, j Node].
        """
        nodes = [list_nodes[0][member[0]], list_nodes[1][member[0]], list_nodes[2][member[0]],
                 list_nodes[0][member[1]], list_nodes[1][member[1]], list_nodes[2][member[1]]]
        self.Members.append(nodes)
        self.ConvertToPrint()

    def AddSupports(self, list_nodes: List[List[float]], support: List[Union[int,bool]] | np.ndarray) -> None:
        """
        Adds supports to be displeyed in the 3D rendering.

        :param list_nodes: list. List of all nodes in the frame. shape: (3, # nodes)
        :param support: list. Support nodes and what DOF are supported. shape: [i Node, DX, DY, DZ, RX, RY, RZ]
        """

        location = [list_nodes[0][int(support[0])], list_nodes[1][int(support[0])], list_nodes[2][int(support[0])]]
        supports = support[1:]
        self.supports.append([location, supports])
        self.ConvertToPrint()

    def AddReleces(self, list_members: List[List[int]], list_nodes: List[List[float]],
                   releces: List[Union[int,bool]]) -> None:
        """
        Adds releces to be displeyed in the 3D rendering.

        :param list_members: list. List of all members in the frame. shape: (2: [i Node, j Node], # nodes)
        :param list_nodes: list. List of all nodes in the frame. shape: (3, # nodes)
        :param releces: list, Releces for a member.
                        shape: [i member, iDX, iDY, iDZ, iRX, iRY, iRZ, jDX, jDY, jDZ, jRX, jRY, jRZ]
        """

        n1 = [list_nodes[0][int(list_members[0][int(releces[0])])],
              list_nodes[1][int(list_members[0][int(releces[0])])],
              list_nodes[2][int(list_members[0][int(releces[0])])]]

        n2 = [list_nodes[0][int(list_members[1][int(releces[0])])],
              list_nodes[1][int(list_members[1][int(releces[0])])],
              list_nodes[2][int(list_members[1][int(releces[0])])]]

        locations = n1 + n2
        direction = releces[1:]
        self.releces.append([locations, direction, int(releces[0])])
        self.ConvertToPrint()

    def AddNodeLoads(self, list_nodes: List[List[float]], node_loads: List[Union[int,float]]) -> None:
        """
        Adds node loads to be displeyed in the 3D rendering.

        :param list_nodes: list. List of all nodes in the frame. shape: (3, # nodes)
        :param node_loads: list. Node load. shape: [i Node, PX, PY, PZ, MX, MY, MZ]
        """

        locations = [list_nodes[0][int(node_loads[0])],
                     list_nodes[1][int(node_loads[0])],
                     list_nodes[2][int(node_loads[0])]]
        load = node_loads[1:]
        self.point_loads.append([locations, load])

    def AddMemberPointLoads(self, list_members: List[List[int]], list_nodes: List[List[float]],
                            member_point_loads: List[Union[int,float]]) -> None:
        """
        Adds member point loads to be displeyed in the 3D rendering.

        :param list_members: list. List of all members in the frame. shape: (2: [i Node, j Node], # nodes)
        :param list_nodes: list. List of all nodes in the frame. shape: (3, # nodes)
        :param member_point_loads: list. Local member point load. shape: [i Member, x, PX, PY, PZ, MX, MY, MZ]
        """
        location = get_loc_on_member(int(member_point_loads[0]), member_point_loads[1], list_members, list_nodes)
        Trans = get_member_t(int(member_point_loads[0]), list_members, list_nodes)[:3, :3]
        P_global = Trans.T @ np.array([member_point_loads[2],member_point_loads[3],member_point_loads[4]]) @ Trans
        M_global = Trans.T @ np.array([member_point_loads[5],member_point_loads[6],member_point_loads[7]]) @ Trans
        self.point_loads.append([location, np.concatenate((P_global,M_global))])

    def AddMemberDistLoads(self, list_members: List[List[int]], list_nodes: List[List[float]],
                           member_dist_loads: List[Union[int,bool]]) -> None:
        """
        Adds member distributed loads to be displeyed in the 3D rendering.

        :param list_members: list. List of all members in the frame. shape: (2: [i Node, j Node], # nodes)
        :param list_nodes: list. List of all nodes in the frame. shape: (3, # nodes)
        :param member_dist_loads: list. Local member distributed load.
                                  shape: [i Member, x1, x2, wx1, wx2, wy1, wy2, wz1, wz2]
        """

        location1 = get_loc_on_member(member_dist_loads[0], member_dist_loads[1], list_members, list_nodes)
        location2 = get_loc_on_member(member_dist_loads[0], member_dist_loads[2], list_members, list_nodes)
        location = location1.tolist() + location2.tolist()
        Trans = get_member_t(member_dist_loads[0], list_members, list_nodes)[:3, :3]
        P1_global = Trans.T @ np.array([member_dist_loads[3], member_dist_loads[5], member_dist_loads[7]]) @ Trans
        P2_global = Trans.T @ np.array([member_dist_loads[4], member_dist_loads[6], member_dist_loads[8]]) @ Trans
        load = [P1_global[0],P2_global[0],P1_global[1],P2_global[1],P1_global[2],P2_global[2]]
        self.dist_loads.append([location, load])

    def ConvertToPrint(self):
        """
        Converts all arrays of things to be displayed into the individual parts that can be printed.
        """

        self.PrintNodes: np.ndarray = np.empty((0, 3))
        self.PrintLines: np.ndarray = np.empty((0, 2, 3))
        self.PrintSurfaceTri: np.ndarray = np.empty((0, 3, 3))
        self.PrintSolidTri: np.ndarray = np.empty((0, 3, 3))
        self.PrintText: list = [[], np.empty((0, 6))]

        # nodes
        for i in range(len(self.Nodes)):
            self.PrintNodes = np.vstack((self.PrintNodes, self.Nodes[i]))

        # releces
        self._convert_releces()

        # members
        for i in range(len(self.Members)):
            self.PrintLines = np.vstack((self.PrintLines,
                                         [[[self.Members[i][0], self.Members[i][1], self.Members[i][2]],
                                          [self.Members[i][3], self.Members[i][4], self.Members[i][5]]]]))

        # supports
        self._convert_supports()

        # point loads
        self._convert_point_loads(self.window.LookDir)

        # distributed loads
        self._convert_dist_loads(self.window.LookDir)

        # text
        for i in range(len(self.text)):
            self.PrintText[0].append(self.text[i][0])
            self.PrintText[1] = np.vstack((self.PrintText[1],[[self.text[i][1],self.text[i][2],self.text[i][3],
                                                                 self.text[i][4],self.text[i][5],self.text[i][6]]]))

    def _convert_supports(self):
        s = self.scale[0]

        X_circle = [s / 4 * np.cos(2 * np.pi * i / self.res[0]) for i in range(self.res[0])]
        Y_circle = [s / 4 * np.sin(2 * np.pi * i / self.res[0]) for i in range(self.res[0])]

        for support in self.supports:
            x = support[0][0]
            y = support[0][1]
            z = support[0][2]

            if all(support[1]): # fully fixed
                # solid horizontal rectangle at top
                n1 = [x - s, y - s, z]
                n2 = [x - s, y + s, z]
                n3 = [x + s, y + s, z]
                n4 = [x + s, y - s, z]
                tri = [[n1, n2, n3], [n1, n4, n3]]
                for i in range(len(tri)): self.PrintSurfaceTri = np.vstack((self.PrintSurfaceTri, [tri[i]]))

            else:
                if support[1][2]:  # vertical support

                    if support[1][0]:  # x horizontal support
                        if support[1][4]:  # y rotation support
                            # x-axis rectangle
                            n1 = [x - s, y, z]
                            n2 = [x - s, y, z - s]
                            n3 = [x + s, y, z - s]
                            n4 = [x + s, y, z]
                            self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n3, n4]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n4, n1]]))

                        else:
                            # x-axis triangle
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

                        if support[1][4]:  # y rotation support
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
                        if support[1][3]:  # x rotation support
                            # y-axis rectangle
                            n1 = [x, y - s, z]
                            n2 = [x, y - s, z - s]
                            n3 = [x, y + s, z - s]
                            n4 = [x, y + s, z]
                            self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n3, n4]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n4, n1]]))
                        else:
                            # y-axis triangle
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

                        if support[1][3]:  # X rotation support
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

                    if support[1][5]:  # Z rotation support
                        # horizontal rectangle under
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

                        if support[1][5]:  # Z rotation support
                            # horizontal half rectangle
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

                        if support[1][5]:  # z rotation support
                            # horizontal half rectangle
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

                    if support[1][3]:  # x rotation support
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

                    if support[1][4]:  # Y rotation support
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

                    if support[1][5] and not (support[1][0] and support[1][1]):  # z rotation support
                        # horizontal rectangle
                        n1 = [x - s, y - s, z]
                        n2 = [x - s, y + s, z]
                        n3 = [x + s, y + s, z]
                        n4 = [x + s, y - s, z]
                        self.PrintLines = np.vstack((self.PrintLines, [[n1, n2]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n2, n3]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n3, n4]]))
                        self.PrintLines = np.vstack((self.PrintLines, [[n4, n1]]))

    def _convert_releces(self) -> None:

        for i in range(len(self.releces)):
            release = np.array(self.releces[i][1])

            n1 = np.array(self.releces[i][0][0:3])
            n2 = np.array(self.releces[i][0][3:6])
            v = n2 - n1
            v /= 10
            r1 = release[0:5]
            r2 = release[6:11]

            if any(r1):
                p1 = n1 + v
                self.PrintNodes = np.vstack((self.PrintNodes, [p1]))

                # noinspection PyTypeChecker
                self.Members[self.releces[i][2]][0] = p1[0]
                # noinspection PyTypeChecker
                self.Members[self.releces[i][2]][1] = p1[1]
                # noinspection PyTypeChecker
                self.Members[self.releces[i][2]][2] = p1[2]

            if any(r2):
                p2 = n2 + v
                self.PrintNodes = np.vstack((self.PrintNodes, [p2]))

                # noinspection PyTypeChecker
                self.Members[self.releces[i][2]][3] = p2[0]
                # noinspection PyTypeChecker
                self.Members[self.releces[i][2]][4] = p2[1]
                # noinspection PyTypeChecker
                self.Members[self.releces[i][2]][5] = p2[2]

    def _convert_point_loads(self, look_dir: np.ndarray) -> None:
        largest: float = 0

        for i in range(len(self.point_loads)):
            for j in range(6):
                if abs(self.point_loads[i][1][j]) > largest:
                    largest = self.point_loads[i][1][j]

        if largest != 0:
            scaler: float  = self.scale[1] / largest

            X_circle = [scaler / 2 * np.cos(2 * np.pi * i / self.res[0]) for i in range(self.res[0])]
            Y_circle = [scaler / 2 * np.sin(2 * np.pi * i / self.res[0]) for i in range(self.res[0])]

            for i in range(len(self.point_loads)):
                for j in range(6):
                    if not np.isclose(self.point_loads[i][1][j],0):
                        a=1
                        p1 = np.array(self.point_loads[i][0])

                        # point load
                        if j <= 2:
                            v = np.zeros(3)
                            v[j] = self.point_loads[i][1][j]
                            p2 = p1 + v  * scaler

                            p3 = p1 + v * scaler / 10

                            unit_v = v / np.linalg.norm(v)

                            side_v = np.cross(unit_v, look_dir)
                            side_v /= np.linalg.norm(side_v)

                            side_v *= scaler / 10

                            p4 = p3 + side_v
                            p3 -= side_v

                            self.PrintLines = np.vstack((self.PrintLines, [[p1, p2]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[p1, p3]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[p1, p4]]))

                            # todo add text label

                        # moment load
                        else:
                            v = np.zeros(3)
                            v[j-3] = self.point_loads[i][1][j] * scaler
                            unit_v = v / np.linalg.norm(v)

                            center = p1 + unit_v * self.scale[1] / 4

                            if unit_v[0] == 1 or unit_v[0] == -1:
                                X = [center[0]] * int(self.res[0] * 0.75)
                                Y = [center[1] + X_circle[i] for i in range(int(self.res[0] * 0.75))]
                                Z = [center[2] + Y_circle[i] for i in range(int(self.res[0] * 0.75))]

                            elif unit_v[1] == 1 or unit_v[1] == -1:
                                X = [center[0] + X_circle[i] for i in range(int(self.res[0] * 0.75))]
                                Y = [center[1]] * int(self.res[0] * 0.75)
                                Z = [center[2] + Y_circle[i] for i in range(int(self.res[0] * 0.75))]

                            elif unit_v[2] == 1 or unit_v[2] == -1:
                                X = [center[0] + X_circle[i] for i in range(int(self.res[0] * 0.75))]
                                Y = [center[1] + Y_circle[i] for i in range(int(self.res[0] * 0.75))]
                                Z = [center[2]] * int(self.res[0] * 0.75)

                            else:
                                raise Exception('Error. Moment load printing invalid moment.')

                            for k in range(int(self.res[0] * 0.75) - 1):
                                self.PrintLines = np.vstack((self.PrintLines, [[[X[k], Y[k], Z[k]],
                                                                                [X[k + 1], Y[k + 1], Z[k + 1]]]]))

                            if unit_v[0] == -1 or unit_v[1] == -1 or unit_v[2] == -1:
                                n1 = np.array([X[0], Y[0], Z[0]])
                                n2 = np.array([X[1], Y[1], Z[1]])
                            else:
                                n1 = np.array([X[len(X) - 1], Y[len(X) - 1], Z[len(X) - 1]])
                                n2 = np.array([X[len(X) - 2], Y[len(X) - 2], Z[len(X) - 2]])

                            v_curve = n2 - n1
                            unit_v_curve = v_curve / np.linalg.norm(v_curve)

                            side_v = np.cross(unit_v,unit_v_curve)

                            n3 = n1 + unit_v_curve * scaler / 10
                            side_v /= np.linalg.norm(side_v)

                            side_v *= scaler / 10

                            n4 = n3 + side_v
                            n3 -= side_v

                            self.PrintLines = np.vstack((self.PrintLines, [[n1, n3]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[n1, n4]]))

                            # todo add text label

    def _convert_dist_loads(self, look_dir: np.ndarray) -> None:

        largest: float = 0

        for i in range(len(self.dist_loads)):
            for j in range(6):
                if abs(self.dist_loads[i][1][j]) > largest:
                    largest = self.dist_loads[i][1][j]

        if largest != 0:
            scaler: float = self.scale[1] / largest

            for i in range(len(self.dist_loads)):
                for j in range(3):

                    if not np.isclose(self.dist_loads[i][1][j*2], 0) and not np.isclose(self.dist_loads[i][1][j*2+1], 0):
                        p1 = np.array([self.dist_loads[i][0][0],self.dist_loads[i][0][1],self.dist_loads[i][0][2]])
                        p2 = np.array([self.dist_loads[i][0][3],self.dist_loads[i][0][4],self.dist_loads[i][0][5]])
                        v = p2 - p1
                        len_v = np.linalg.norm(v)

                        num_inter = math.ceil(len_v / self.scale[2])
                        dist_between = len_v / num_inter

                        v_gap = v / len_v * dist_between

                        v_force_1 = self.dist_loads[i][1][j*2]  * scaler
                        v_force_2 = self.dist_loads[i][1][j*2+1]  * scaler

                        p3 = np.copy(p1)
                        p4 = np.copy(p2)

                        p3[j] += v_force_1
                        p4[j] += v_force_2

                        v_top = p4 - p3
                        len_v_top = np.linalg.norm(v_top)

                        dif_v_force = v_force_2 - v_force_1

                        change_v_force = dif_v_force / num_inter

                        v_gap_top = np.copy(v_gap)
                        v_gap_top[j] += change_v_force

                        self.PrintLines = np.vstack((self.PrintLines, [[p3, p4]]))

                        for k in range(num_inter+1):
                            p_bot = p1 + v_gap * k
                            p_top = p3 + v_gap_top * k
                            v_line = p_top - p_bot

                            # arow
                            unit_v_line = v_line / np.linalg.norm(v_line)

                            p_arow1 = p_bot + unit_v_line * scaler / 10

                            side_v_line = np.cross(unit_v_line, look_dir)
                            side_v_line /= np.linalg.norm(side_v_line)

                            side_v_line *= scaler / 10

                            p_arow2 = p_arow1 + side_v_line
                            p_arow1 -= side_v_line

                            self.PrintLines = np.vstack((self.PrintLines, [[p_bot, p_top]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[p_bot, p_arow1]]))
                            self.PrintLines = np.vstack((self.PrintLines, [[p_bot, p_arow2]]))

                        # todo add text label

    # todo add converter to convert everything to the print arrays
    # todo rework all the adders
    # todo popup to select what should be printed

def get_loc_on_member(i_member: int, x, list_members: List[List[int]], list_nodes: List[List[float]]):
    """
    Finds the global coords of where x lands on a member.

    :param i_member: int, index of member.
    :param x: float. location along member.
    :param list_members: list. List of all members in the frame. shape: (2: [i Node, j Node], # nodes)
    :param list_nodes: list. List of all nodes in the frame. shape: (3, # nodes)
    :return: Location on member.
    """

    n1 = np.array([list_nodes[0][int(list_members[0][int(i_member)])],
                   list_nodes[1][int(list_members[0][int(i_member)])],
                   list_nodes[2][int(list_members[0][int(i_member)])]])

    n2 = np.array([list_nodes[0][int(list_members[1][int(i_member)])],
                   list_nodes[1][int(list_members[1][int(i_member)])],
                   list_nodes[2][int(list_members[1][int(i_member)])]])

    vector = n2 - n1
    unit_vector = vector / np.linalg.norm(vector)

    return  n1 + unit_vector * x

def get_member_t(i_member: int, list_members: List[List[int]], list_nodes: List[List[float]]):
    """
    Builds an array of the transformation matrices for the member.

    :param i_member: int, index of member.
    :param list_members: list. List of all members in the frame. shape: (2: [i Node, j Node], # nodes)
    :param list_nodes: list. List of all nodes in the frame. shape: (3, # nodes)
    :return: transformation matrix.
    """

    # Get the global coordinates for the two ends
    n1 = np.array([list_nodes[0][int(list_members[0][int(i_member)])],
                   list_nodes[1][int(list_members[0][int(i_member)])],
                   list_nodes[2][int(list_members[0][int(i_member)])]])

    n2 = np.array([list_nodes[0][int(list_members[1][int(i_member)])],
                   list_nodes[1][int(list_members[1][int(i_member)])],
                   list_nodes[2][int(list_members[1][int(i_member)])]])

    i = list_members[0][int(i_member)]
    j = list_members[1][int(i_member)]
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