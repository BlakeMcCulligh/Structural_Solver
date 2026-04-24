
class WrongNumberOfBoundsError(Exception):
    pass

class InvalidMemberTypeError(Exception):
    pass

class CrossSection:
    def __init__(self):
        """
        Initializes the cross-section object.
        """
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

    def setOptBoundsForFrames(self, minBounds: list, maxBounds: list):
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


    def add2DTruss(self, A: float, E: float):
        """
        Assignes the values for a 2D truss
        :param A: Area
        :param E: Elastic Modulus
        """
        self.A = A
        self.E = E

    def add2DFrame(self, A: float, E: float, I: float):
        """
        Assignes the values for a 2D frame
        :param A: Area
        :param E: Elastic Modulus
        :param I: Stifness Modulus
        """
        self.A = A
        self.E = E
        self.I_main = I

    def add3DTruss(self, A: float, E: float):
        """
        Assignes the values for a 3D truss
        :param A: Area
        :param E: Elastic Modulus
        """
        self.A = A
        self.E = E

    def add3DFrame(self, A: float, E: float, G: float, I_main: float, I_weak: float, J: float):
        """
        Assignes the values for a 3D frame
        :param A: Area
        :param E: Elastic Modulus
        :param G: shear Modulus
        :param I_main: Stifness Modulus for the strong axis
        :param I_weak: Stifness Modulus for the weak axis
        :param J: Stifness Modulus for torsion
        """
        self.A = A
        self.E = E
        self.G = G
        self.I_main = I_main
        self.I_weak = I_weak
        self.J = J

    def optAdd3DTruss(self, E: float, minBounds: float, maxBounds: float):
        """
        Assignes the values for a 3D truss to be optimized
        :param E: Elastic Modulus
        :param minBounds: minimum area
        :param maxBounds: maximum area
        """
        self.E = E
        self.minBounds = [minBounds]
        self.maxBounds = [maxBounds]

    def optAdd3DFrame(self, E: float, G: float, memberType: str, minBounds: list, maxBounds: list):
        """
        Assignes the values for a 3D frame to be optimized
        :param E: Elastic Modulus
        :param G: shear Modulus
        :param memberType: Type of member
        :param minBounds: minimum bounds on the optimization variables
        :param maxBounds: maximum bounds on the optimization variables
        """
        self.E = E
        self.G = G
        self.Type = memberType

        self.setOptBoundsForFrames(minBounds, maxBounds)

    def optAdd2DTruss(self, E: float, minBounds: float, maxBounds: float):
        """
        Assignes the values for a 2D truss to be optimized
        :param E: Elastic Modulus
        :param minBounds: minimum area
        :param maxBounds: maximum area
        """
        self.E = E
        self.minBounds = [minBounds]
        self.maxBounds = [maxBounds]

    def optAdd2DFrame(self, E: float, memberType: str, minBounds: list, maxBounds: list):
        """
        Assignes the values for a 2D frame to be optimized
        :param E: Elastic Modulus
        :param memberType: Type of member
        :param minBounds: minimum bounds on the optimization variables
        :param maxBounds: maximum bounds on the optimization variables
        """
        self.E = E
        self.Type = memberType

        self.setOptBoundsForFrames(minBounds, maxBounds)