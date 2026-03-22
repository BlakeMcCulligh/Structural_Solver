
class WrongNumberOfCordsGivenError(Exception):
    pass

class WrongNumberOfDegressOfFredomGivenError(Exception):
    pass

class Node:
    def __init__(self, isTruss: bool, is3D: bool, cords: list):
        self.isTruss = isTruss
        self.is3D = is3D

        self.cords = None

        if self.is3D:
            if len(cords) == 3:
                self.cords = cords
            else:
                raise WrongNumberOfCordsGivenError("Number of cords given must be 3 for 3D structures")
        else:
            if len(cords) == 2:
                self.cords = cords
            else:
                raise WrongNumberOfCordsGivenError("Number of cords given must be 2 for 2D structures")

        self.support = None
        self.load = None

        if self.isTruss:
            if self.is3D:
                self.support = [False, False, False]
                self.load = [0, 0, 0]
            else:
                self.support = [False, False]
                self.load = [0, 0]
        else:
            if self.is3D:
                self.support = [False, False, False, False, False, False]
                self.load = [0, 0, 0, 0, 0, 0]
            else:
                self.support = [False, False, False]
                self.load = [0, 0, 0]

    def addSupport(self, support: list):
        if self.isTruss:
            if self.is3D:
                if len(support) == 3:
                    self.support = support
                else:
                    raise WrongNumberOfDegressOfFredomGivenError("Number of degres given must be 3 for a 3D truss")
            else:
                if len(support) == 2:
                    self.support = support
                else:
                    raise WrongNumberOfDegressOfFredomGivenError("Number of degres given must be 2 for a 2D truss")
        else:
            if self.is3D:
                if len(support) == 6:
                    self.support = support
                else:
                    raise WrongNumberOfDegressOfFredomGivenError("Number of degres given must be 6 for a 3D frame")
            else:
                if len(support) == 3:
                    self.support = support
                else:
                    raise WrongNumberOfDegressOfFredomGivenError("Number of degres given must be 3 for a 2D frame")

    def addLoad(self, load: list):
        if self.isTruss:
            if self.is3D:
                if len(load) == 3:
                    self.load = load
                else:
                    raise WrongNumberOfDegressOfFredomGivenError("Number of degres given must be 3 for a 3D truss")
            else:
                if len(load) == 2:
                    self.load = load
                else:
                    raise WrongNumberOfDegressOfFredomGivenError("Number of degres given must be 2 for a 2D truss")
        else:
            if self.is3D:
                if len(load) == 6:
                    self.load = load
                else:
                    raise WrongNumberOfDegressOfFredomGivenError("Number of degres given must be 6 for a 3D frame")
            else:
                if len(load) == 3:
                    self.load = load
                else:
                    raise WrongNumberOfDegressOfFredomGivenError(f"Number of degres given must be 3 for a 2D frame. {len(load)} were given.")





