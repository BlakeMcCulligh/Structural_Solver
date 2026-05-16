"""
Holds the object that stores all the 3D frame's data
and handles distributing data when a new part of the frame is added or edited.
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

from frame_3D_gui.display import get_loc_on_member, get_member_t


class Data:
    """
    Frame Object holding all the data needed to define the frame.
    """

    def __init__(self):
        """
        Frame Object Constructor.
        """

        self.Nodes = [[], [], []]
        self.Materials = [[], [], [], [], []]
        self.Members = [[], [], [], [], [], [], [], []]
        self.Supports = [[], [], [], [], [], [], []]
        self.Releases = [[], [], [], [], [], [], [], [], [], [], [], [], []]
        self.NodeLoad = [[], [], [], [], [], [], [], []]
        self.MemberPointLoad = [[], [], [], [], [], [], [], [], []]
        self.MemberDistLoad = [[], [], [], [], [], [], [], [], [], []]

    def AddNodes(self, Window, Nodes, AddToTables: bool, AddToDisplay: bool):
        """
        Adds nodes to the frame.

        :param Window: Object storing the main window.
        :param Nodes: List of node coordinates to be added to the frame.
        :param AddToTables: If the nodes are to be added to the input table.
        :param AddToDisplay: If the nodes are to be added to the 3D display.
        """

        if isinstance(Nodes[0], float):
            Nodes[0] = [Nodes[0]]
            Nodes[1] = [Nodes[1]]
            Nodes[2] = [Nodes[2]]

        self.Nodes[0] = self.Nodes[0] + Nodes[0]
        self.Nodes[1] = self.Nodes[1] + Nodes[1]
        self.Nodes[2] = self.Nodes[2] + Nodes[2]

        if AddToTables:
            num_row = len(Window.Tables[0].get_children())
            for i in range(len(Nodes[0])):
                Window.Tables[0].insert('', 'end', values=[str(num_row + i), Nodes[0][i], Nodes[1][i], Nodes[2][i]])

        if AddToDisplay:
            Nodes = np.array(Nodes)
            for i in range(len(Nodes[0])):
                #Window.AddPrintNode(Nodes[:, i])
                Window.DisplayData.add_node(Nodes[:, i])

        _reset_solutions(Window)

    def EditNode(self, Window, NewNode, Index, EditTable, EditDisplay, RowID):
        """
        Edits the coordinates of a node in the frame.

        :param Window: Object storing the main window.
        :param NewNode: New coordinates for the node.
        :param Index: Index of the node to edit.
        :param EditTable: If the input Tables are to be edited.
        :param EditDisplay: If the 3D display is to be edited.
        :param RowID: the row within the input table that the edited node is in.
        """

        for i in range(3): self.Nodes[i][Index] = NewNode[i]

        if EditTable:
            col_ids = Window.Tables[0]["columns"]
            Window.Tables[0].set(RowID, col_ids[0], Index)
            for i in range(3): Window.Tables[0].set(RowID, col_ids[i + 1], NewNode[i])

        if EditDisplay:
            Window.DisplayData.Nodes[Index] = NewNode

            # updating members
            Window.DisplayData.Members = []
            Members = np.array(self.Members)
            for i in range(len(self.Members[0])):
                Window.DisplayData.add_member(self.Nodes, Members[[0, 1], i].astype(int))

            # updating supports
            Window.DisplayData.supports = []
            Supports = np.array(self.Supports)
            for i in range(len(self.Supports[0])):
                Window.DisplayData.add_supports(self.Nodes, Supports[:, i])

            Window.DisplayData.convert_to_print()
            Window.update_canvas()

        _reset_solutions(Window)

    def AddMaterials(self, Window, Material, AddToTables: bool, AddToDisplay: bool):
        """
        Adds materials to the frame.

        :param Window: Object storing the main window.
        :param Material: list of materials to be Added to the frame.
        :param AddToTables: If the input tables are to be edited.
        :param AddToDisplay: If the 3D display is to be edited.
        """

        for i in range(5):
            if isinstance(Material[i], float): Material[i] = [Material[i]]
            self.Materials[i] = self.Materials[i] + Material[i]

        if AddToTables:
            num_row = len(Window.Tables[1].get_children())
            for i in range(len(Material[0])):
                Window.Tables[1].insert('', 'end', values=[str(num_row + i), Material[0][i], Material[1][i],
                                                           Material[2][i], Material[3][i], Material[4][i]])

        if AddToDisplay:
            pass #TODO

        _reset_solutions(Window)

    def EditMaterials(self, Window, NewMaterials, Index, EditTable, EditDisplay, RowID):
        """
        Edits a material in the frame.

        :param Window:  Object storing the main window.
        :param NewMaterials: New property values for the material.
        :param Index: Index of the material to be edited.
        :param EditTable: If the input tables are to be edited.
        :param EditDisplay: If the 3D display is to be edited.
        :param RowID: the row within the input table that the edited material is in.
        """

        for i in range(5): self.Materials[i][Index] = NewMaterials[i]

        if EditTable:
            col_ids = Window.Tables[1]["columns"]
            Window.Tables[1].set(RowID, col_ids[0], Index)
            for i in range(5): Window.Tables[1].set(RowID, col_ids[i + 1], NewMaterials[i])

        if EditDisplay:
            pass # todo

        _reset_solutions(Window)

    def AddMembers(self, Window, Members, AddToTables: bool, AddToDisplay: bool):
        """
        Adds members to the frame.

        :param Window:  Object storing the main window.
        :param Members: List of members to be Added to the frame.
        :param AddToTables: If the input tables are to be edited.
        :param AddToDisplay: If the 3D display is to be edited.
        """

        for i in range(8):
            if isinstance(Members[i], float)  or isinstance(Members[i], bool): Members[i] = [Members[i]]
            self.Members[i] = self.Members[i] + Members[i]

        if AddToTables:
            num_row = len(Window.Tables[2].get_children())
            for i in range(len(Members[0])):
                Window.Tables[2].insert('', 'end', values=[str(num_row + i), Members[0][i], Members[1][i],
                                                           Members[2][i], Members[3][i], Members[4][i], Members[5][i],
                                                           Members[6][i], Members[7][i]])

        if AddToDisplay:
            Members = np.array(Members)
            for i in range(len(Members[0])):
                #Window.AddPrintLine(Members[[0, 1], i])
                Window.DisplayData.add_member(self.Nodes, Members[[0, 1], i].astype(int))

        _reset_solutions(Window)

    def EditMembers(self, Window, NewMembers, Index, EditTable, EditDisplay, RowID):
        """
        Edits a member in the frame.

        :param Window: Object storing the main window.
        :param NewMembers: New property values for the member.
        :param Index: Index of the member to be edited.
        :param EditTable: If the input tables are to be Edited.
        :param EditDisplay: If the 3D display is to be Edited.
        :param RowID: The input table row containing the member being edited.
        """

        for i in range(8):
            self.Nodes[i][Index] = NewMembers[i]

        if EditTable:
            col_ids = Window.Tables[2]["columns"]
            Window.Tables[2].set(RowID, col_ids[0], Index)
            for i in range(8): Window.Tables[2].set(RowID, col_ids[i + 1], NewMembers[i])

        if EditDisplay:
            Window.DisplayData.Member[Index] = NewMembers[0, 1]
            Window.DisplayData.convert_to_print()
            Window.update_canvas()

        _reset_solutions(Window)

    def AddSupports(self, Window, Supports, AddToTables: bool, AddToDisplay: bool):
        """
        Adds supports to the frame.

        :param Window: Object storing the main window.
        :param Supports: List of Supports to be Added to the frame.
        :param AddToTables: If the input tables are to be edited.
        :param AddToDisplay: If the 3D display is to be edited.
        """

        for i in range(7):
            if isinstance(Supports[i], float) or isinstance(Supports[i], bool): Supports[i] = [Supports[i]]
            self.Supports[i] = self.Supports[i] + Supports[i]

        if AddToTables:
            num_row = len(Window.Tables[3].get_children())
            for i in range(len(Supports[0])):
                Window.Tables[3].insert('', 'end', values=[str(num_row + i), Supports[0][i], Supports[1][i],
                                                           Supports[2][i], Supports[3][i], Supports[4][i],
                                                           Supports[5][i], Supports[6][i]])

        if AddToDisplay:
            for i in range(len(Supports[0])):
                sup = []
                for j in range(7):
                    sup.append(Supports[j][i])
                Window.DisplayData.AddSupports(self.Nodes, sup)
            Window.UpdateCanves()

        _reset_solutions(Window)

    def EditSupports(self, Window, NewSupport, Index, EditTable, EditDisplay, RowID):
        """
        Edits a supports in the frame.

        :param Window: Object storing the main window.
        :param NewSupport: The new values for the support being edited.
        :param Index: The Index of the support being edited.
        :param EditTable: If the input tables are to be edited.
        :param EditDisplay: If the 3D display is to be edited.
        :param RowID: The input table row containing the support being edited.
        """

        for i in range(7): self.Supports[i][Index] = NewSupport[i]

        if EditTable:
            col_ids = Window.Tables[3]["columns"]
            Window.Tables[3].set(RowID, col_ids[0], Index)
            for i in range(7): Window.Tables[3].set(RowID, col_ids[i + 1], NewSupport[i])

        if EditDisplay:
            sup = []
            for j in range(7): sup.append(self.Supports[j][Index])
            sup = np.array(sup)
            location = [self.Nodes[0][int(sup[0])], self.Nodes[1][int(sup[0])], self.Nodes[2][int(sup[0])]]
            supports = sup[1:7]
            Window.DisplayData.supports[Index] = [location, supports]
            Window.DisplayData.convert_to_print()
            Window.update_canvas()

        _reset_solutions(Window)

    def AddReleases(self, Window, Releases, AddToTables: bool, AddToDisplay: bool):
        """
        Adds releases to the frame.

        :param Window: Object storing the main window.
        :param Releases: List of releases to be Added to the frame.
        :param AddToTables: If the input tables are to be edited.
        :param AddToDisplay: If the 3D display is to be edited.
        """

        for i in range(13):
            if isinstance(Releases[i], float) or isinstance(Releases[i], bool): Releases[i] = [Releases[i]]
            self.Releases[i] = self.Releases[i] + Releases[i]

        if AddToTables:
            num_row = len(Window.Tables[4].get_children())
            for i in range(len(Releases[0])):
                Window.Tables[4].insert('', 'end', values=[str(num_row + i), Releases[0][i], Releases[1][i],
                                                           Releases[2][i], Releases[3][i], Releases[4][i],
                                                           Releases[5][i], Releases[6][i], Releases[7][i],
                                                           Releases[8][i], Releases[9][i], Releases[10][i],
                                                           Releases[11][i], Releases[12][i]])

        if AddToDisplay:
            for i in range(len(Releases[0])):
                release = []
                for j in range(13): release.append(self.Releases[j][i])
                Window.DisplayData.AddReleces([self.Members[0],self.Members[1]], self.Nodes, release)

            Window.DisplayData.ConvertToPrint()
            Window.UpdateCanves()

        _reset_solutions(Window)

    def EditReleases(self, Window, NewReleases, Index, EditTable, EditDisplay, RowID):
        """
        Edits a releases in the frame.

        :param Window: Object storing the main window.
        :param NewReleases: The New values for the release being edited.
        :param Index: Index of the release being edited.
        :param EditTable: If the input tables are to be edited.
        :param EditDisplay: If the 3D display is to be edited.
        :param RowID: The input table row containing the release being edited.
        """

        for i in range(12): self.Materials[i][Index] = NewReleases[i]

        if EditTable:
            col_ids = Window.Tables[4]["columns"]
            Window.Tables[4].set(RowID, col_ids[0], Index)
            for i in range(12): Window.Tables[4].set(RowID, col_ids[i + 1], NewReleases[i])

        if EditDisplay:
            n1 = [self.Nodes[0][self.Members[0][NewReleases[0]]],
                  self.Nodes[1][self.Members[0][NewReleases[0]]],
                  self.Nodes[2][self.Members[0][NewReleases[0]]]]
            n2 = [self.Nodes[0][self.Members[1][NewReleases[0]]],
                  self.Nodes[1][self.Members[1][NewReleases[0]]],
                  self.Nodes[2][self.Members[1][NewReleases[0]]]]
            locations = n1 + n2
            direction = NewReleases[1:]
            Window.DisplayData.releces[Index] = [locations, direction]
            Window.DisplayData.ConvertToPrint()
            Window.UpdateCanves()

        _reset_solutions(Window)

    def AddNodeLoads(self, Window, NodeLoad, AddToTables: bool, AddToDisplay: bool):
        """
        Adds Node point loads to the frame.

        :param Window: Object storing the main window.
        :param NodeLoad: List of node point loads to be added to the frame.
        :param AddToTables: If the input tables are to be edited.
        :param AddToDisplay: If the 3D display is to be edited.
        """

        for i in range(8):
            if isinstance(NodeLoad[i], float): NodeLoad[i] = [NodeLoad[i]]
            self.NodeLoad[i] = self.NodeLoad[i] + NodeLoad[i]

        if AddToTables:
            num_row = len(Window.Tables[5].get_children())
            for i in range(len(NodeLoad[0])):
                Window.Tables[5].insert('', 'end', values=[str(num_row + i), NodeLoad[0][i], NodeLoad[1][i],
                                                           NodeLoad[2][i], NodeLoad[3][i], NodeLoad[4][i],
                                                           NodeLoad[5][i], NodeLoad[6][i], NodeLoad[7][i]])

        if AddToDisplay:
            for i in range(len(NodeLoad[0])):
                load = []
                for j in range(8): load.append(self.NodeLoad[j][i])
                Window.DisplayData.AddNodeLoads(self.Nodes, load)

            Window.DisplayData.ConvertToPrint()
            Window.UpdateCanves()

        _reset_solutions(Window)

    def EditNodeLoads(self, Window, NewNodeLoad, Index, EditTable, EditDisplay, RowID):
        """
        Edits a node point loads in the frame.

        :param Window: Object storing the main window.
        :param NewNodeLoad: The new values for the node load being edited.
        :param Index: Index of the Node load being edited.
        :param EditTable: If the input tables are to be edited.
        :param EditDisplay: If the 3D display is to be edited.
        :param RowID: The input table row containing the node load being edited.
        """

        for i in range(8): self.Materials[i][Index] = NewNodeLoad[i]

        if EditTable:
            col_ids = Window.Tables[5]["columns"]
            Window.Tables[5].set(RowID, col_ids[0], Index)
            for i in range(8): Window.Tables[5].set(RowID, col_ids[i + 1], NewNodeLoad[i])

        if EditDisplay:
            location = [self.Nodes[0][Index],self.Nodes[1][Index],self.Nodes[2][Index]]
            load = NewNodeLoad[1:]
            Window.DisplayData.point_loads[Index] = [location, load]
            Window.DisplayData.ConvertToPrint()
            Window.UpdateCanves()

        _reset_solutions(Window)

    def AddMemberPointLoads(self, Window, MemberPointLoad, AddToTables: bool, AddToDisplay: bool):
        """
        Adds member point loads to the frame.

        :param Window: Object storing the main window.
        :param MemberPointLoad: List of member point loads to be Added to the frame.
        :param AddToTables: If the input Tables are to be edited.
        :param AddToDisplay: If the 3D display is to be edited.
        """

        for i in range(9):
           if isinstance(MemberPointLoad[i], float): MemberPointLoad[i] = [MemberPointLoad[i]]
           self.MemberPointLoad[i] = self.MemberPointLoad[i] + MemberPointLoad[i]

        if AddToTables:
           num_row = len(Window.Tables[6].get_children())
           for i in range(len(MemberPointLoad[0])):
               Window.Tables[6].insert('', 'end', values=[str(num_row + i), MemberPointLoad[0][i], MemberPointLoad[1][i],
                                                          MemberPointLoad[2][i], MemberPointLoad[3][i],
                                                          MemberPointLoad[4][i], MemberPointLoad[5][i],
                                                          MemberPointLoad[6][i], MemberPointLoad[7][i],
                                                          MemberPointLoad[8][i]])

        if AddToDisplay:
            for i in range(len(MemberPointLoad[0])):
               load = []
               for j in range(9): load.append(self.MemberPointLoad[j][i])
               Window.DisplayData.AddMemberPointLoads([self.Members[0],self.Members[1]], self.Nodes, load)

            Window.DisplayData.ConvertToPrint()
            Window.UpdateCanves()

        _reset_solutions(Window)

    def EditMemberPointLoad(self, Window, NewMemberPointLoad, Index, EditTable, EditDisplay, RowID):
        """
        Edits a member point loads in the frame.

        :param Window: Object storing the main window.
        :param NewMemberPointLoad: The new values for the member point load being edited.
        :param Index: Index of the member point load being edited.
        :param EditTable: If the input tables are to be edited.
        :param EditDisplay: If the 3D display is to be edited.
        :param RowID: The Input table row containing the member point load being edited.
        """

        for i in range(9): self.Materials[i][Index] = NewMemberPointLoad[i]

        if EditTable:
            col_ids = Window.Tables[6]["columns"]
            Window.Tables[6].set(RowID, col_ids[0], Index)
            for i in range(9): Window.Tables[6].set(RowID, col_ids[i + 1], NewMemberPointLoad[i])

        if EditDisplay:
            location = get_loc_on_member(int(NewMemberPointLoad[0]), int(NewMemberPointLoad[1]), self.Members, self.Nodes)
            Trans = get_member_t(int(NewMemberPointLoad[0]), self.Members, self.Nodes)[:3, :3]
            P_global = Trans.T @ np.array([NewMemberPointLoad[2], NewMemberPointLoad[3], NewMemberPointLoad[4]]) @ Trans
            M_global = Trans.T @ np.array([NewMemberPointLoad[5], NewMemberPointLoad[6], NewMemberPointLoad[7]]) @ Trans
            Window.DisplayData.point_loads = [location, np.concatenate((P_global,M_global))]
            Window.DisplayData.ConvertToPrint()
            Window.UpdateCanves()

        _reset_solutions(Window)

    def AddMemberDistLoads(self, Window, MemberDistLoad, AddToTables: bool, AddToDisplay: bool):
        """
        Adds member distributed loads to the frame.

        :param Window:  Object storing the main window.
        :param MemberDistLoad: List of member distributed loads to be Added to the frame.
        :param AddToTables: If the input tables are to be edited.
        :param AddToDisplay: If the 3D display is to be edited.
        """

        for i in range(10):
            if isinstance(MemberDistLoad[i], float): MemberDistLoad[i] = [MemberDistLoad[i]]
            self.MemberDistLoad[i] = self.MemberDistLoad[i] + MemberDistLoad[i]

        if AddToTables:
            num_row = len(Window.Tables[7].get_children())
            for i in range(len(MemberDistLoad[0])):
                Window.Tables[7].insert('', 'end', values=[str(num_row + i), MemberDistLoad[0][i], MemberDistLoad[1][i],
                                                           MemberDistLoad[2][i], MemberDistLoad[3][i],
                                                           MemberDistLoad[4][i], MemberDistLoad[5][i],
                                                           MemberDistLoad[6][i], MemberDistLoad[7][i],
                                                           MemberDistLoad[8][i], MemberDistLoad[9][i]])

        if AddToDisplay:
            for i in range(len(MemberDistLoad[0])):
                load = []
                for j in range(10): load.append(self.MemberDistLoad[j][i])
                Window.DisplayData.AddMemberDistLoads([self.Members[0],self.Members[1]], self.Nodes, load)

            Window.DisplayData.ConvertToPrint()
            Window.UpdateCanves()

        _reset_solutions(Window)

    def EditMemberDistLoad(self, Window, NewMemberDistLoad, Index, EditTable, EditDisplay, RowID):
        """
        Edits a member distributed loads in the frame.

        :param Window:  Object storing the main window.
        :param NewMemberDistLoad: The new values for the member distributed load being edited.
        :param Index: Index of the member distributed load being edited.
        :param EditTable: If the input tables are to be edited.
        :param EditDisplay: If the 3D display is to be edited.
        :param RowID: The input table row containing the member distributed load being edited.
        """

        for i in range(10): self.Materials[i][Index] = NewMemberDistLoad[i]

        if EditTable:
            col_ids = Window.Tables[7]["columns"]
            Window.Tables[7].set(RowID, col_ids[0], Index)
            for i in range(10): Window.Tables[7].set(RowID, col_ids[i + 1], NewMemberDistLoad[i])

        if EditDisplay:
            location1 = get_loc_on_member(NewMemberDistLoad[0], NewMemberDistLoad[1], self.Members, self.Nodes)
            location2 = get_loc_on_member(NewMemberDistLoad[0], NewMemberDistLoad[2], self.Members, self.Nodes)
            location = location1 + location2
            Trans = get_member_t(NewMemberDistLoad[0], self.Members, self.Nodes)[:3, :3]
            P1_global = Trans.T @ np.array([NewMemberDistLoad[3], NewMemberDistLoad[5], NewMemberDistLoad[7]]) @ Trans
            P2_global = Trans.T @ np.array([NewMemberDistLoad[4], NewMemberDistLoad[6], NewMemberDistLoad[8]]) @ Trans
            load = [P1_global[0], P2_global[0], P1_global[1], P2_global[1], P1_global[2], P2_global[2]]
            Window.DisplayData.dist_loads[Index] = [location, load]
            Window.DisplayData.ConvertToPrint()
            Window.UpdateCanves()

        _reset_solutions(Window)

def _reset_solutions(window):
    """
    Sets all solutions to not solved, when anything is changed within the frame.

    :param window: Object storing the main window.
    """

    window.Frame = None
    window.Results = None
    window.OptimizationResults = None