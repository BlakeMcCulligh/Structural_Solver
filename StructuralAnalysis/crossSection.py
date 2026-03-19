
class CrossSection:
    def __init__(self):
        self.isTruss = False
        self.is3D = False

        self.memberType = None

        self.A = None
        self.E = None
        self.G = None
        self.I_main = None
        self.I_weak = None
        self.J = None

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

    def optAdd3DFrame(self, E: float, G: float, memberType: str):
        self.E = E
        self.G = G
        self.memberType = memberType

    def optAdd2DTruss(self, E: float):
        self.E = E

    def optAdd2DFrame(self, E: float, memberType: str):
        self.E = E
        self.memberType = memberType