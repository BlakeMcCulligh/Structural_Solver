
class CrossSection:
    def __init__(self):
        self.A = None # Area
        self.Ag = None # Gross Area

        self.Ax = None # effective cross-sectional area for tention
        self.Asy = None # efective shear area in the y direction
        self.Asz = None # efective shear area in the z direction

        self.Iz = None
        self.Iy = None
        self.J = None

        self.Cw = None
        self.numSymmetry = None

        self.rx = None
        self.ry = None
        self.xo = None
        self.yo = None