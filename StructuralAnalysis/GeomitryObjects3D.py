
class Element:
    def __init__(self, N1, N2):
        self.N1 = N1
        self.N2 = N2
        self.L = ((N2[0]-N1[0])**2 + (N2[1]-N1[1])**2 + (N2[2]-N1[2])**2)**0.5