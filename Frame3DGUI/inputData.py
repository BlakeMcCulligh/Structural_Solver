import numpy as np


class data:
    def __init__(self):

        self.nodes = [[],[],[]]
        self.materials = [[],[],[],[],[]]
        self.members = [[],[],[],[],[],[],[],[]]
        self.supports = [[],[],[],[],[],[],[]]
        self.releases = [[],[],[],[],[],[],[],[],[],[],[],[],[]]
        self.nodeLoad = [[],[],[],[],[],[],[],[]]
        self.memberPointLoad = [[],[],[],[],[],[],[],[],[]]
        self.memberDistLoad = [[],[],[],[],[],[],[],[],[],[]]

    def addnodes(self, window, nodes, addToTables: bool, addToDisplay: bool):

        if isinstance(nodes[0], float):
            nodes[0] = [nodes[0]]
            nodes[1] = [nodes[1]]
            nodes[2] = [nodes[2]]

        self.nodes[0] = self.nodes[0] + nodes[0]
        self.nodes[1] = self.nodes[1] + nodes[1]
        self.nodes[2] = self.nodes[2] + nodes[2]

        if addToTables:
            numRow = len(window.tables[0].get_children())
            for i in range(len(nodes[0])):
                window.tables[0].insert('', 'end', values=[str(numRow + i), nodes[0][i], nodes[1][i], nodes[2][i]])


        if addToDisplay:
            nodes = np.array(nodes)
            for i in range(len(nodes[0])):
                window.addPrintNode(nodes[:, i])

        resetSolutions(window)

    def editnode(self, window, newNode, index,  editTable, editDisplay, row_id):
        for i in range(3): self.nodes[i][index] = newNode[i]

        if editTable:
            col_ids = window.tables[0]["columns"]
            window.tables[0].set(row_id, col_ids[0], index)
            for i in range(3): window.tables[0].set(row_id, col_ids[i+1], newNode[i])

        if editDisplay:
            window.printNodes[index] = newNode
            window.updateCanves()

        resetSolutions(window)

    def addmaterials(self, window, material, addToTables: bool, addToDisplay: bool):
        for i in range(5):
            if isinstance(material[i], float): material[i] = [material[i]]
            self.materials[i] = self.materials[i] + material[i]

        if addToTables:
            numRow = len(window.tables[1].get_children())
            for i in range(len(material[0])):
                window.tables[1].insert('', 'end', values=[str(numRow + i), material[0][i], material[1][i], material[2][i], material[3][i], material[4][i]])

        if addToDisplay:
            pass #TODO

        resetSolutions(window)

    def editmaterials(self, window, newMaterials, index, editTable, editDisplay, row_id):
        for i in range(5): self.materials[i][index] = newMaterials[i]

        if editTable:
            col_ids = window.tables[1]["columns"]
            window.tables[1].set(row_id, col_ids[0], index)
            for i in range(5): window.tables[1].set(row_id, col_ids[i + 1], newMaterials[i])

        if editDisplay:
            pass # todo

        resetSolutions(window)

    def addmembers(self, window, members, addToTables: bool, addToDisplay: bool):

        for i in range(8):
            if isinstance(members[i], float)  or isinstance(members[i], bool): members[i] = [members[i]]
            self.members[i] = self.members[i] + members[i]

        if addToTables:
            numRow = len(window.tables[2].get_children())
            for i in range(len(members[0])):
                window.tables[2].insert('', 'end', values=[str(numRow + i), members[0][i], members[1][i],
                                                           members[2][i], members[3][i],members[4][i], members[5][i],
                                                           members[6][i], members[7][i]])

        if addToDisplay:
            members = np.array(members)
            for i in range(len(members[0])):
                window.addPrintLine(members[[0, 1], i])

        resetSolutions(window)

    def editmembers(self, window, newMembers, index, editTable, editDisplay, row_id):
        for i in range(8):
            self.nodes[i][index] = newMembers[i]

        if editTable:
            col_ids = window.tables[2]["columns"]
            window.tables[2].set(row_id, col_ids[0], index)
            for i in range(8): window.tables[2].set(row_id, col_ids[i + 1], newMembers[i])

        if editDisplay:
            window.printLine[index] = newMembers[0, 1]
            window.updateCanves()

        resetSolutions(window)

    def addsupports(self, window, supports, addToTables: bool, addToDisplay: bool):
        for i in range(7):
            if isinstance(supports[i], float) or isinstance(supports[i], bool): supports[i] = [supports[i]]
            self.supports[i] = self.supports[i] + supports[i]

        if addToTables:
            numRow = len(window.tables[3].get_children())
            for i in range(len(supports[0])):
                window.tables[3].insert('', 'end', values=[str(numRow + i), supports[0][i], supports[1][i],
                                                           supports[2][i], supports[3][i], supports[4][i],
                                                           supports[5][i], supports[6][i]])

        if addToDisplay:
            pass #TODO

        resetSolutions(window)

    def editsupports(self, window, newSupport, index, editTable, editDisplay, row_id):
        for i in range(7): self.materials[i][index] = newSupport[i]

        if editTable:
            col_ids = window.tables[3]["columns"]
            window.tables[3].set(row_id, col_ids[0], index)
            for i in range(7): window.tables[3].set(row_id, col_ids[i + 1], newSupport[i])

        if editDisplay:
            pass # todo

        resetSolutions(window)

    def addreleases(self, window, releases, addToTables: bool, addToDisplay: bool):
        for i in range(13):
            if isinstance(releases[i], float) or isinstance(releases[i], bool): releases[i] = [releases[i]]
            self.releases[i] = self.releases[i] + releases[i]

        if addToTables:
            numRow = len(window.tables[4].get_children())
            for i in range(len(releases[0])):
                window.tables[4].insert('', 'end', values=[str(numRow + i), releases[0][i], releases[1][i],
                                                           releases[2][i], releases[3][i],releases[4][i],
                                                           releases[5][i], releases[6][i], releases[7][i],
                                                           releases[8][i],releases[9][i], releases[10][i],
                                                           releases[11][i], releases[12][i]])

        if addToDisplay:
            pass #TODO

        resetSolutions(window)

    def editreleases(self, window, newReleases, index, editTable, editDisplay, row_id):
        for i in range(12): self.materials[i][index] = newReleases[i]

        if editTable:
            col_ids = window.tables[4]["columns"]
            window.tables[4].set(row_id, col_ids[0], index)
            for i in range(12): window.tables[4].set(row_id, col_ids[i + 1], newReleases[i])

        if editDisplay:
            pass # todo

        resetSolutions(window)

    def addnodeLoads(self, window, nodeLoad, addToTables: bool, addToDisplay: bool):
        for i in range(8):
            if isinstance(nodeLoad[i], float): nodeLoad[i] = [nodeLoad[i]]
            self.nodeLoad[i] = self.nodeLoad[i] + nodeLoad[i]

        if addToTables:
            numRow = len(window.tables[5].get_children())
            for i in range(len(nodeLoad[0])):
                window.tables[5].insert('', 'end', values=[str(numRow + i), nodeLoad[0][i], nodeLoad[1][i],
                                                           nodeLoad[2][i], nodeLoad[3][i], nodeLoad[4][i],
                                                           nodeLoad[5][i], nodeLoad[6][i], nodeLoad[7][i]])

        if addToDisplay:
            pass #TODO

        resetSolutions(window)

    def editnodeLoads(self, window, newNodeLoad, index, editTable, editDisplay, row_id):
        for i in range(8): self.materials[i][index] = newNodeLoad[i]

        if editTable:
            col_ids = window.tables[5]["columns"]
            window.tables[5].set(row_id, col_ids[0], index)
            for i in range(8): window.tables[5].set(row_id, col_ids[i + 1], newNodeLoad[i])

        if editDisplay:
            pass # todo

        resetSolutions(window)

    def addmemberPointLoads(self, window, memberPointLoad, addToTables: bool, addToDisplay: bool):
       for i in range(9):
           if isinstance(memberPointLoad[i], float): memberPointLoad[i] = [memberPointLoad[i]]
           self.memberPointLoad[i] = self.memberPointLoad[i] + memberPointLoad[i]

       if addToTables:
           numRow = len(window.tables[6].get_children())
           for i in range(len(memberPointLoad[0])):
               window.tables[6].insert('', 'end', values=[str(numRow + i), memberPointLoad[0][i], memberPointLoad[1][i],
                                                          memberPointLoad[2][i], memberPointLoad[3][i],
                                                          memberPointLoad[4][i], memberPointLoad[5][i],
                                                          memberPointLoad[6][i], memberPointLoad[7][i],
                                                          memberPointLoad[8][i]])

       if addToDisplay:
           pass #TODO

       resetSolutions(window)

    def editmemberPointLoad(self, window, newMemberPointLoad, index, editTable, editDisplay, row_id):
        for i in range(9): self.materials[i][index] = newMemberPointLoad[i]

        if editTable:
            col_ids = window.tables[6]["columns"]
            window.tables[6].set(row_id, col_ids[0], index)
            for i in range(9): window.tables[6].set(row_id, col_ids[i + 1], newMemberPointLoad[i])

        if editDisplay:
            pass # todo

        resetSolutions(window)

    def addmemberDistLoads(self, window, memberDistLoad, addToTables: bool, addToDisplay: bool):
        for i in range(10):
            if isinstance(memberDistLoad[i], float): memberDistLoad[i] = [memberDistLoad[i]]
            self.memberDistLoad[i] = self.memberDistLoad[i] + memberDistLoad[i]

        if addToTables:
            numRow = len(window.tables[7].get_children())
            for i in range(len(memberDistLoad[0])):
                window.tables[7].insert('', 'end', values=[str(numRow + i), memberDistLoad[0][i], memberDistLoad[1][i],
                                                           memberDistLoad[2][i], memberDistLoad[3][i],
                                                           memberDistLoad[4][i], memberDistLoad[5][i],
                                                           memberDistLoad[6][i], memberDistLoad[7][i],
                                                           memberDistLoad[8][i],memberDistLoad[9][i]])

        if addToDisplay:
            pass #TODO

        resetSolutions(window)

    def editmemberDistLoad(self, window, newMemberDistLoad, index, editTable, editDisplay, row_id):
        for i in range(10): self.materials[i][index] = newMemberDistLoad[i]

        if editTable:
            col_ids = window.tables[7]["columns"]
            window.tables[7].set(row_id, col_ids[0], index)
            for i in range(10): window.tables[7].set(row_id, col_ids[i + 1], newMemberDistLoad[i])

        if editDisplay:
            pass # todo

        resetSolutions(window)

def resetSolutions(window):
    window.Frame = None
    window.results = None
    window.Optimization_Results = None