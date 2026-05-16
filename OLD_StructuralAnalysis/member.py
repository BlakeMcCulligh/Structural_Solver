
class WrongNumberOfDegressOfFredomGivenError(Exception):
    pass

class Member:
    def __init__(self, nodes: list, crossSections: list, is3D: bool, isTruss: bool, nodeIndexs: list, crossSectionIndex: int):
        """
        Creation of a member object.
        :param nodes: the list of Node objectsin the structure.
        :param crossSections: the list of crossSection objects in the structure.
        :param is3D: is the member a part of a 3D structure.
        :param isTruss: is the member a part of a truss.
        :param nodeIndexs: the indexes of the PrintNodes to be the start and end of the member.
        :param crossSectionIndex: the i of the crossSection of the memebr.
        """
        self.isTruss = isTruss
        self.is3D = is3D

        self.nodeIndexs = nodeIndexs

        self.nodes = [nodes[self.nodeIndexs[0]], nodes[self.nodeIndexs[1]]]

        self.crossSectionIndex = crossSectionIndex
        self.crossSection = crossSections[self.crossSectionIndex]

        if self.isTruss:
            self.relece = None
        else:
            if self.is3D:
                self.relece = [False] * 12
            else:
                self.relece = [False] * 6


        # Loads
        self.pointLoadsMagnatude = []
        self.pointLoadsLocation = []

        self.uniformlyDistributedLoads = []

        self.distributedLoadsMagnatude = []
        self.distributedLoadsLocation = []

    def addRelece(self, relece: list):
        """
        Updates the member's releases to the paramiter. Only for Frames.
        :param relece: the new member's releases.
        """
        if not self.isTruss:
            if self.is3D:
                if len(relece) == 12:
                    self.relece = relece
                else:
                    raise WrongNumberOfDegressOfFredomGivenError("Number of degres given must be 12 for a 3D frame Releces")
            else:
                if len(relece) == 6:
                    self.relece = relece
                else:
                    raise WrongNumberOfDegressOfFredomGivenError("Number of degres given must be 6 for a 2D frame Releces")
        else:
            print("releases can not be set for trusses")

    def addPointLoad(self, magnatude: list, location: float):
        if not self.isTruss:
            if self.is3D:
                if len(magnatude) == 6:
                    self.pointLoadsMagnatude.append(magnatude)
                    self.pointLoadsLocation.append(location)
                else:
                    raise WrongNumberOfDegressOfFredomGivenError("Number of degres given must be 6 for a 3D frame point loads")
            else:
                if len(magnatude) == 3:
                    self.pointLoadsMagnatude.append(magnatude)
                    self.pointLoadsLocation.append(location)
                else:
                    raise WrongNumberOfDegressOfFredomGivenError("Number of degres given must be 3 for a 3D frame point loads")
        else:
            if self.is3D:
                if len(magnatude) == 3:
                    self.pointLoadsMagnatude.append(magnatude)
                    self.pointLoadsLocation.append(location)
                else:
                    raise WrongNumberOfDegressOfFredomGivenError("Number of degres given must be 3 for a 3D truss point loads")
            else:
                if len(magnatude) == 2:
                    self.pointLoadsMagnatude.append(magnatude)
                    self.pointLoadsLocation.append(location)
                else:
                    raise WrongNumberOfDegressOfFredomGivenError("Number of degres given must be 2 for a 2D truss point loads")

    def addUniformlyDistributedLoad(self, magnatude: list):
        if self.is3D:
            if len(magnatude) == 3:
                self.uniformlyDistributedLoads.append(magnatude)
            else:
                raise WrongNumberOfDegressOfFredomGivenError("Number of degres given must be 3 for a 3D structre uniformly distributed loads")
        else:
            if len(magnatude) == 2:
                self.uniformlyDistributedLoads.append(magnatude)
            else:
                raise WrongNumberOfDegressOfFredomGivenError("Number of degres given must be 2 for a 2D structre uniformly distributed loads")

    def addDistributedLoad(self, magnatude: list, location: list):
        if self.is3D:
            if len(magnatude) == 6:
                if len(location) == 2:
                    self.distributedLoadsMagnatude.append(magnatude)
                    self.distributedLoadsLocation.append(location)
                else:
                    raise WrongNumberOfDegressOfFredomGivenError("Number of distances given must be 2 for distributed loads")
            else:
                raise WrongNumberOfDegressOfFredomGivenError("Number of degres given must be 6 for a 3D structre distributed loads")
        else:
            if len(magnatude) == 4:
                if len(location) == 2:
                    self.distributedLoadsMagnatude.append(magnatude)
                    self.distributedLoadsLocation.append(location)
                else:
                    raise WrongNumberOfDegressOfFredomGivenError("Number of distances given must be 2 for distributed loads")
            else:
                raise WrongNumberOfDegressOfFredomGivenError("Number of degres given must be 4 for a 2D structre distributed loads")

    def getEndReactions(self):
        #todo
        pass

