"""
Used for cheacking the acercy of the 3D frame solver.
"""

from Pynite import FEModel3D

from frame_3D_solver.main import Frame3D

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

# TODO add comments

class TestCase:
    def __init__(self):
        self.nodes = []
        self.members = []
        self.materials = []
        self.sections = []
        self.supports = []
        self.releases = []
        self.num_casses = 0
        self.node_loads = []
        self.member_point_loads = []
        self.member_dist_loads = []

    def AddNode(self, x,y,z):
        self.nodes.append([x,y,z])

    def AddMaterial(self, E, G, nu, rho):
        self.materials.append([E,G,nu,rho])

    def AddSection(self, A, Iy, Iz, J):
        self.sections.append([A,Iy,Iz,J])

    def AddMember(self, n1, n2, materoal, section):
        self.members.append([n1, n2, materoal, section])

    def AddSupport(self, node, fx, fy, fz, mx, my, mz):
        self.supports.append([node, fx, fy, fz, mx, my, mz])

    def AddRelease(self, member, ifx, ify, ifz, imx, imy, imz, jfx, jfy, jfz, jmx, jmy, jmz):
        self.releases.append([member, ifx, ify, ifz, imx, imy, imz, jfx, jfy, jfz, jmx, jmy, jmz])

    def AddNodeLoad(self, node, fx, fy, fz, mx, my, mz, case):
        self.node_loads.append([node, fx, fy, fz, mx, my, mz, case])

    def AddMemberPointLoad(self, member, x, fx, fy, fz, mx, my, mz, case):
        self.member_point_loads.append([member, x, fx, fy, fz, mx, my, mz, case])

    def AddMemberDistLoad(self, member, x1, x2, wx1, wx2, wy1, wy2, wz1, wz2, case):
        self.member_dist_loads.append([member, x1, x2, wx1, wx2, wy1, wy2, wz1, wz2, case])

    def SetNumCasses(self, num_casses):
        self.num_casses = num_casses

    def Check(self):

        results_pynite = self._get_results_PyNite()
        print("Pynite: ", results_pynite)

        results_structsolver = self._get_results_struct_solver()
        print("StructSolver: ", results_structsolver)

        # TODO Add method to compair the resoults and flage any miss matches

    def _get_results_PyNite(self):
        model = FEModel3D()

        for i, node in enumerate(self.nodes):
            model.add_node(str(i), node[0], node[1], node[2])

        for i, material in enumerate(self.materials):
            model.add_material(str(i), material[0], material[1], material[2], material[3])

        for i, section in enumerate(self.sections):
            model.add_section(str(i), section[0], section[1], section[2], section[3])

        for i, member in enumerate(self.members):
            model.add_member(str(i), str(member[0]), str(member[1]), str(member[2]), str(member[3]))

        for i, support in enumerate(self.supports):
            model.def_support(str(support[0]), bool(support[1]), bool(support[2]), bool(support[3]), bool(support[4]),
                              bool(support[5]), bool(support[6]))

        for i, release in enumerate(self.releases):
            model.def_releases(str(release[0]), bool(release[1]), bool(release[2]), bool(release[3]), bool(release[4]),
                               bool(release[5]), bool(release[6]), bool(release[7]), bool(release[8]), bool(release[9]),
                               bool(release[10]), bool(release[11]), bool(release[12]))

        directions = ['FX', 'FY', 'FZ', 'MX', 'MY', 'MZ']
        for i, nl in enumerate(self.node_loads):
            for j in range(6):
                if nl[j + 1] != 0: model.add_node_load(str(nl[0]), direction=directions[j], P=nl[j + 1], case=str(nl[7]))

        for i, mpl in enumerate(self.member_point_loads):
            for j in range(6):
                if mpl[j + 2] != 0:
                    model.add_member_pt_load(str(mpl[0]), x=mpl[1], direction=directions[j], P=mpl[j + 2], case=str(mpl[8]))

        for i, mdl in enumerate(self.member_dist_loads):
            if mdl[3] != 0 or mdl[4] != 0:
                model.add_member_dist_load(str(mdl[0]), 'FX', mdl[3], mdl[4], mdl[1], mdl[2], str(mdl[9]))
            if mdl[5] != 0 or mdl[6] != 0:
                model.add_member_dist_load(str(mdl[0]), 'FY', mdl[5], mdl[6], mdl[1], mdl[2], str(mdl[9]))
            if mdl[7] != 0 or mdl[8] != 0:
                model.add_member_dist_load(str(mdl[0]), 'FZ', mdl[7], mdl[8], mdl[1], mdl[2], str(mdl[9]))

        for i in range(self.num_casses):
            model.add_load_combo(str(i), {str(i):1}, combo_tags=None)

        model.analyze_linear(log=False)

        node_deflections = []
        for i in range(len(self.nodes)):
            for j in range(self.num_casses):
                node_deflections.append([model.nodes[str(i)].DX[str(j)],model.nodes[str(i)].DY[str(j)],
                                         model.nodes[str(i)].DZ[str(j)],model.nodes[str(i)].RX[str(j)],
                                         model.nodes[str(i)].RY[str(j)],model.nodes[str(i)].RZ[str(j)]])

        return node_deflections

    def _get_results_struct_solver(self):

        model = Frame3D(None)

        for node in self.nodes:
            model.AddNode(node[0], node[1], node[2])

        for material in self.materials:
            model.AddMaterial(material[0], material[1], material[2], material[3])

        for member in self.members:
            model.AddMember(member[0], member[1], member[2], True, self.sections[member[3]][0],
                            self.sections[member[3]][1], self.sections[member[3]][2], self.sections[member[3]][3])

        for support in self.supports:
            model.AddSupport(support[0], support[1], support[2], support[3], support[4], support[5], support[6])

        for release in self.releases:
            model.AddReleases(release[0], release[1], release[2], release[3], release[4], release[5], release[6],
                              release[7], release[8], release[9], release[10], release[11], release[12])

        for nl in self.node_loads:
            model.AddNodeLoad(nl[0],nl[1],nl[2],nl[3],nl[4],nl[5],nl[6],nl[7])

        for mpl in self.member_point_loads:
            model.AddMemberPointLoad(mpl[0], mpl[1], mpl[2], mpl[3], mpl[4], mpl[5], mpl[6], mpl[7], mpl[8])

        for mdl in self.member_dist_loads:
            model.AddMemberDistLoad(mdl[0], mdl[1], mdl[2], mdl[3], mdl[4], mdl[5], mdl[6], mdl[7], mdl[8], mdl[9])

        model.PreAnalysisLinear()
        D, DX, DY, DZ, RX, RY, RZ, weight, reactions, internal_forces = model.AnalysisLinear()

        return D, DX, DY, DZ, RX, RY, RZ, weight, reactions, internal_forces

t = TestCase()
t.AddNode(0,0,0)
t.AddNode(1,0,0)

t.AddMaterial(1,1,1,1)
t.AddSection(1,1,1,1)

t.AddMember(0,1,0,0)

t.AddSupport(0,1,1,1,1,1,1)

t.AddNodeLoad(1,1,1,1,0,0,0,0)
#t.AddMemberPointLoad(0,0.1,1,0,0,0,0,0,0)
#t.AddMemberDistLoad(0,0,1,1,0,0,0,0,0,0)

t.SetNumCasses(1)

t.Check()