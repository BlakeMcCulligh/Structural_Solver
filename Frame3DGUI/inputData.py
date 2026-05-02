
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

    def addnodes(self, nodes):
        self.nodes[0] = self.nodes[0] + nodes[0]
        self.nodes[1] = self.nodes[1] + nodes[1]
        self.nodes[2] = self.nodes[2] + nodes[2]

    def addmaterials(self, material):
        for i in range(5):
            self.materials[i] = self.materials[i] + material[i]

    def addmembers(self, members):
        for i in range(8):
            self.members[i] = self.members[i] + members[i]

    def addsupports(self, supports):
        for i in range(7):
            self.supports[i] = self.supports[i] + supports[i]

    def addreleases(self, releases):
        for i in range(13):
            self.releases[i] = self.releases[i] + releases[i]

    def addnodeLoads(self, nodeLoad):
        for i in range(8):
            self.nodeLoad[i] = self.nodeLoad[i] + nodeLoad[i]

    def addmemberPointLoads(self, memberPointLoad):
       for i in range(9):
           self.memberPointLoad[i] = self.memberPointLoad[i] + memberPointLoad[i]

    def addmemberDistLoads(self, memberDistLoad):
        for i in range(10):
            self.memberDistLoad[i] = self.memberDistLoad[i] + memberDistLoad[i]