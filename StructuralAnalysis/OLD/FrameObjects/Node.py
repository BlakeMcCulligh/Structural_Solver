
class Node:
    def __init__(self):
        self.numDimensions = None
        self.cords = []

    def calcNumDimensions(self):
        self.numDimensions = len(self.cords)
