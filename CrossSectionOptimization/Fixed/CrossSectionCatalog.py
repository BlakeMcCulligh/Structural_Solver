
class Catalog:
    def __init__(self, sections):
        self.sections = sections
        self.n = len(sections)

class CrossSection:
    def __init__(self, A, I, rho, sigma_y):
        self.A = A
        self.I = I
        self.rho = rho
        self.sigma_y = sigma_y
