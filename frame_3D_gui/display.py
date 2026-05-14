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

        self.scale = []

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

    def ConvertToPrint(self):
        """
        Converts all arrays of things to be displayed into the individual parts that can be printeded.
        """

        self.PrintNodes: np.ndarray = np.empty((0, 3))
        self.PrintLines: np.ndarray = np.empty((0, 2, 3))
        self.PrintSurfaceTri: np.ndarray = np.empty((0, 3, 3))
        self.PrintSolidTri: np.ndarray = np.empty((0, 3, 3))

        for i in range(len(self.Nodes)):
            self.PrintNodes = np.vstack((self.PrintNodes, self.Nodes[i]))

        for i in range(len(self.Members)):
            self.PrintLines = np.vstack((self.PrintLines,
                                         [[[self.Members[i][0], self.Members[i][1], self.Members[i][2]],
                                          [self.Members[i][3], self.Members[i][4], self.Members[i][5]]]]))


    # todo add inputs for arrays
    # todo add converter to convert everything to the print arrays
    # todo rework all the adders
    # todo popup to select what should be printed