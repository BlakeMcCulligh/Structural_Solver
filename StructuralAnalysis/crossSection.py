
class WrongNumberOfBoundsError(Exception):
    pass

class InvalidMemberTypeError(Exception):
    pass

class CrossSection:
    def __init__(self):
        self.isTruss = False
        self.is3D = False

        self.Type = None

        self.A = None
        self.E = None
        self.G = None
        self.I_main = None
        self.I_weak = None
        self.J = None

        self.minBounds = None
        self.maxBounds = None

    def add2DTruss(self, A: float, E: float):
        self.A = A
        self.E = E

    def add2DFrame(self, A: float, E: float, I: float):
        self.A = A
        self.E = E
        self.I_main = I

    def add3DTruss(self, A: float, E: float):
        self.A = A
        self.E = E

    def add3DFrame(self, A: float, E: float, G: float, I_main: float, I_weak: float, J: float):
        self.A = A
        self.E = E
        self.G = G
        self.I_main = I_main
        self.I_weak = I_weak
        self.J = J

    def optAdd3DTruss(self, E: float):
        self.E = E

    def optAdd3DFrame(self, E: float, G: float, memberType: str, minBounds: list, maxBounds: list):
        self.E = E
        self.G = G
        self.Type = memberType
        self.minBounds = minBounds
        self.maxBounds = maxBounds

    def optAdd2DTruss(self, E: float):
        self.E = E

    def optAdd2DFrame(self, E: float, memberType: str, minBounds: list, maxBounds: list):
        self.E = E
        self.Type = memberType

        if self.Type == "SquareHSS" or self.Type == "TubeHSS":
            if len(minBounds) == 2:
                self.minBounds = minBounds
            else:
                WrongNumberOfBoundsError(f"Wrong number of bounds. Length of minBounds should be 2 for {self.Type}")
            if len(maxBounds) == 2:
                self.maxBounds = maxBounds
            else:
                WrongNumberOfBoundsError(f"Wrong number of bounds. Length of maxBounds should be 2 for {self.Type}")
        elif self.Type == "RectHSS" or self.Type == "Angle":
            if len(minBounds) == 3:
                self.minBounds = minBounds
            else:
                WrongNumberOfBoundsError(f"Wrong number of bounds. Length of minBounds should be 3 for {self.Type}")
            if len(maxBounds) == 3:
                self.maxBounds = maxBounds
            else:
                WrongNumberOfBoundsError(f"Wrong number of bounds. Length of maxBounds should be 3 for {self.Type}")
        else:
            InvalidMemberTypeError(f"{self.Type} is not a valid member type. Chose from SquareHSS, RectHSS, TubeHSS, or Angle.")