import numpy as np

class structure:
    def __init__(self):
        self.nodes = []
        self.members = []
        self.loads = {}

    def add_node(self, x, y, fixed=(False, False, False)):
        self.nodes.append(Node(x, y, fixed))
        return len(self.nodes) - 1

    def add_member(self, n1, n2):
        self.members.append(Member(n1, n2))

    def add_load(self, node, fx=0.0, fy=0.0, mz=0.0):
        self.loads[node] = np.array([fx, fy, mz])

class Node:
    def __init__(self, x, y, fixed=(False, False, False)):
        self.x = x
        self.y = y
        self.fixed = fixed  # (ux, uy, rz)

class Member:
    def __init__(self, n1, n2):
        self.n1 = n1
        self.n2 = n2

    def length(self, nodes):
        dx = nodes[self.n2].x - nodes[self.n1].x
        dy = nodes[self.n2].y - nodes[self.n1].y
        return np.hypot(dx, dy)

    def angle(self, nodes):
        dx = nodes[self.n2].x - nodes[self.n1].x
        dy = nodes[self.n2].y - nodes[self.n1].y
        return np.arctan2(dy, dx)