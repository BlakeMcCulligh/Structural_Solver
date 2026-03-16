
class WrongNumberOfDegressOfFredomGivenError(Exception):
    pass

class Member:
    def __init__(self, nodes: list, crossSections: list, is3D: bool, isTruss: bool, nodeIndexs: list, crossSectionIndex: int):
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


    def addRelece(self, relece: list):
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
            print("releces can not be set for trusses")


