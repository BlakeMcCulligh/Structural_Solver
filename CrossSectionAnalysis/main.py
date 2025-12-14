import math

from CrossSectionAnalysis.integration import areaIntegral
from Sketch.geometric_primitives.arc import Arc
from Sketch.geometric_primitives.point import Point
from Sketch.geometric_primitives.segment import Segment
from Sketch.main import sketch
from GeomitryHelpers.lineOrArcIntersect import checkIfLineOrArcIntersect
from GeomitryHelpers.pointWithinShape import point_in_shape

def a(y, x):
    return 1

def analysisCrossSection():
    geometry = sketch()

    cloasedShape = getClosedShape(geometry)

    if cloasedShape is not None:
        Area = cloasedShape.AreaInnerRemoved
        print("points: ")
        for p in cloasedShape.Points:
            print("(", p.x, ",", p.y, ")")
        print(Area)

        answer, errorEstimite = areaIntegral(a, cloasedShape.Geomitry)
        print("Area Testing: ", answer)

    else:
        print("No shape found")

def getClosedShape(geometry):
    closedShapes = findClosedShapes(geometry)

    for shapeOuter in closedShapes:
        for shapeInner in closedShapes:
            if shapeOuter != shapeInner:
                if shapeWitinShape(shapeOuter, shapeInner):
                    shapeOuter.ClosedShapesWithin.append(shapeInner)

    for shapes in closedShapes:
        shapes.updateArea()

    largestShape = None
    if len(closedShapes) != 0:
        largestShape = closedShapes[0]
        for shape in closedShapes:
            if shape.AreaInnerRemoved > largestShape.Area:
                largestShape = shape

    return largestShape

class ClosedShape:
    def __init__(self, geometry, lines, arcs, points, ClosedShapesWithin):
        self.Geomitry = geometry
        self.Lines = lines
        self.Arcs = arcs
        self.Points = points
        self.ClosedShapesWithin = ClosedShapesWithin

        self.Area = self.getArea()
        self.AreaInnerRemoved = None

    def getArea(self):
        area = 0

        g1 = self.Geomitry[0]
        g2 = self.Geomitry[1]

        flip1 = False
        if g1.p1.x == g2.p1.x and g1.p1.y == g2.p1.y:
            flip1 = True
        elif g1.p2.x == g2.p1.x and g1.p2.y == g2.p1.y:
            flip1 = False
        elif g1.p1.x == g2.p2.x and g1.p1.y == g2.p2.y:
            flip1 = True
        elif g1.p2.x == g2.p2.x and g1.p2.y == g2.p2.y:
            flip1 = False
        else:
            print("ERROR")

        lastPoint = [0,0]
        firstG = True
        for g in self.Geomitry:
            if firstG:
                firstG = False
                if flip1:
                    area += 0.5 * (g.p2.x * g.p1.y - g.p1.x * g.p2.y)
                    lastPoint = [g.p1.x, g.p1.y]
                    if isinstance(g, Arc):
                        aArc = self.areaArc(g)
                        addAreaSign = self.findArcAreaSign(g)
                        area += addAreaSign*aArc
                else:
                    area += 0.5 * (g.p1.x * g.p2.y - g.p2.x * g.p1.y)
                    lastPoint = [g.p2.x, g.p2.y]
                    if isinstance(g, Arc):
                        aArc = self.areaArc(g)
                        addAreaSign = self.findArcAreaSign(g)
                        area += addAreaSign * aArc
            else:
                if g.p1.x == lastPoint[0] and g.p1.y == lastPoint[1]:
                    lastPoint = [g.p2.x, g.p2.y]
                    area += 0.5*(g.p1.x*g.p2.y-g.p2.x*g.p1.y)
                    if isinstance(g, Arc):
                        aArc = self.areaArc(g)
                        addAreaSign = self.findArcAreaSign(g)
                        area += addAreaSign * aArc
                else:
                    lastPoint = [g.p1.x, g.p1.y]
                    area += 0.5*(g.p2.x*g.p1.y-g.p1.x*g.p2.y)
                    if isinstance(g, Arc):
                        aArc = self.areaArc(g)
                        addAreaSign = self.findArcAreaSign(g)
                        area += addAreaSign * aArc

        # TODO area with shapes with an arc works

        return abs(area)

    def updateArea(self):
        self.AreaInnerRemoved = self.Area
        for innerShapes in self.ClosedShapesWithin:
            self.AreaInnerRemoved -= innerShapes.Area

    def areaArc(self, arc):
        AreaCircleSection = 0.5 * arc.Radius ** 2 * arc.angle
        AreaTriangle = 0.5 * arc.Radius ** 2 * math.sin(arc.angle)
        return AreaCircleSection - AreaTriangle

    def findArcAreaSign(self, arc):

        midlePoint = arc.middle_point()

        m = (arc.p2.x + arc.p1.x) / (arc.p2.y - arc.p1.y)
        b = arc.p1.y - m * arc.p1.x

        y = m * midlePoint.x + b

        if midlePoint.y > y:
            return 1
        elif midlePoint.y < y:
            return -1
        else:
            print("ERROR")
            return None



def findClosedShapes(geometry):
    loops = []

    edges = geometry.segments + geometry.arcs

    for i in range(len(edges)):
        loopFound = False
        p = []
        e = []

        # adding the first line to their arrays
        p.append(edges[i].p1)
        p.append(edges[i].p2)
        e.append(i)
        loopFound, p, e = loop(i, p, e, loopFound, edges)

        if loopFound:
            loopPoints = p.copy()
            loopGeomitry = []
            for a in e:
                loopGeomitry.append(edges[a])

            if confermLoopNotAlreadyFound(loopGeomitry, loops):
                loopLines, loopArcs = seperateLinesAndArcs(loopGeomitry)
                l = ClosedShape(loopGeomitry, loopLines, loopArcs, loopPoints, [])
                loops.append(l)

    return loops

def loop(i, p, e, loopFound, edges):

    # if a loop is found exit out of the loop function
    if p[0].x == p[len(p)-1].x and p[0].y == p[len(p)-1].y and len(e)>2 or loopFound:
        loopFound = True
    else:
        for j in range(len(edges)):
            # looking for matching points
            if edges[j].p1.x == p[len(p)-1].x and edges[j].p1.y == p[len(p)-1].y:

                # making sure it is not its own point
                notUsed = True
                for m in range(len(e)):
                    if j == e[m]:
                        notUsed = False

                if notUsed:

                    # add point on the other side of the newly found connecting line
                    p.append(edges[j].p2)

                    e.append(j)
                    loopFound, p, e = loop(i, p, e, loopFound, edges)

            if edges[j].p2.x == p[len(p) - 1].x and edges[j].p2.y == p[len(p) - 1].y:

                # making sure it is not its own point
                notUsed = True
                for m in range(len(e)):
                    if j == e[m]:
                        notUsed = False

                if notUsed:
                    # add point on the other side of the newly found connecting line
                    p.append(edges[j].p1)

                    e.append(j)
                    loopFound, p, e = loop(i, p, e, loopFound, edges)

    return loopFound, p, e

def confermLoopNotAlreadyFound(loopGeomitry, loops):
    for loop in loops:
        numSame = 0
        for g in loopGeomitry:
            if isinstance(g, Segment):
                for l in loop.Lines:
                    if linesEqual(l, g):
                        numSame += 1
            else:
                for a in loop.Arcs:
                    if arcsEqual(a, g):
                        numSame += 1
        if numSame == len(loopGeomitry):
            return False
    return True

def linesEqual(l1, l2):
    p1Same = l1.p1.x == l2.p1.x and l1.p1.y == l2.p1.y
    p2Same = l1.p2.x == l2.p2.x and l1.p2.y == l2.p2.y
    return p1Same and p2Same

def arcsEqual(a1, a2):
    p1Same = a1.p1.x == a2.p1.x and a1.p1.y == a2.p1.y
    p2Same = a1.p2.x == a2.p2.x and a1.p2.y == a2.p2.y
    cSame = a1.Center.x == a2.Center.x and a1.Center.y == a2.Center.y
    rSame = a1.radius == a2.radius
    return p1Same and p2Same and cSame and rSame

def seperateLinesAndArcs(loopGeomitry):
    Lines = []
    Arcs = []
    for g in loopGeomitry:
        if isinstance(g, Segment):
            Lines.append(g)
        else:
            Arcs.append(g)

    return Lines, Arcs

def shapeWitinShape(shape1, shape2):
    numPWithin = 0
    for p in shape2.Points:
        if  point_in_shape(p, shape1):
            numPWithin += 1
    if numPWithin >= len(shape2.Points)-1:
        return False

    e1s = shape1.Lines + shape1.Arcs
    e2s = shape2.Lines + shape2.Arcs
    for e1 in e1s:
        for e2 in e2s:
            if checkIfLineOrArcIntersect(e1, e2):
                return False

    return True

analysisCrossSection()
