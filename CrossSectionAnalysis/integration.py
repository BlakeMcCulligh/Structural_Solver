import math
import scipy.integrate as integrate

from Sketch.geometric_primitives.arc import Arc
from Sketch.geometric_primitives.common import v2v_angle_ccw_underPi
from Sketch.geometric_primitives.point import Point
from Sketch.geometric_primitives.segment import Segment
from Sketch.geometric_primitives.vector import Vector

def convertEdgesToEquations(edges):
    bounds = []
    for e in edges:
        if isinstance(e, Segment):
            # line equation
            if round((e.p2.x - e.p1.x)) != 0:
                m = (e.p2.y - e.p1.y) / (e.p2.x - e.p1.x)
                b = e.p1.y - m * e.p2.x

                if m == 0:
                    equation = lambda x, b=b: b
                else:
                    equation = lambda x, m=m, b=b: m * x + b

                if e.p1.x < e.p2.x:
                    bounds.append([equation, e.p1, e.p2, "Line", b, m])
                else:
                    bounds.append([equation, e.p2, e.p1, "Line", b, m])
            else:
                bounds.append([None, e.p1, e.p2, "Line"])

        elif isinstance(e, Arc):
            # arc equation

            equationPlus = lambda x: e.Center.y + (e.Radius**2 - (x - e.Center.x) ** 2)**0.5
            equationMinus = lambda x: e.Center.y - (e.Radius**2 - (x - e.Center.x) ** 2)**0.5

            y1 = equationPlus(e.p1.x)
            y2 = equationPlus(e.p2.x)

            p1InPlus = False
            p2InPlus = False
            if y1 == e.p1.y:
                p1InPlus = True
            if y2 == e.p2.y:
                p2InPlus = True

            if p1InPlus and p2InPlus:
                # add equationPlus
                if e.p1.x < e.p2.x:
                    bounds.append([equationPlus, e.p1, e.p2, "Arc", e.Center.y, e.Radius, e.Center.x, "plus"])
                else:
                    bounds.append([equationPlus, e.p2, e.p1, "Arc", e.Center.y, e.Radius, e.Center.x, "plus"])
            elif not p1InPlus and not p2InPlus:
                # add equationMinus
                if e.p1.x < e.p2.x:
                    bounds.append([equationMinus, e.p1, e.p2, "Arc", e.Center.y, e.Radius, e.Center.x, "Minus"])
                else:
                    bounds.append([equationMinus, e.p2, e.p1, "Arc", e.Center.y, e.Radius, e.Center.x, "Minus"])
            else:

                # add both
                xcPlus = e.Center.x + e.Radius
                xcMinus = e.Center.x - e.Radius
                yc = e.Center.y

                swapPoints = [Point(xcMinus, yc), Point(xcPlus, yc)]

                positiveXVecter = Vector.from_two_points(e.Center, Point(e.Center.x+1, e.Center.y))

                endP1Vector = Vector.from_two_points(e.Center, e.p1)
                endP2Vector = Vector.from_two_points(e.Center, e.p2)

                rangeAngles = [v2v_angle_ccw_underPi(positiveXVecter, endP1Vector), v2v_angle_ccw_underPi(positiveXVecter, endP2Vector)]

                swapP1Vector = Vector.from_two_points(e.Center, swapPoints[0])
                swapP2Vector = Vector.from_two_points(e.Center, swapPoints[1])

                swapAngles = [v2v_angle_ccw_underPi(positiveXVecter, swapP1Vector), v2v_angle_ccw_underPi(positiveXVecter, swapP2Vector)]

                swap1InRange = rangeAngles[0] < swapAngles[0] < rangeAngles[1]
                swap2InRange = rangeAngles[0] < swapAngles[1] < rangeAngles[1]

                if e.p1.x < e.p2.x:
                    p1 = e.p1
                    p2 = e.p2
                else:
                    p1 = e.p2
                    p2 = e.p1

                if swap1InRange and swap2InRange:
                    if p1.y < e.Center.y:
                        bounds.append([equationMinus, p1, swapPoints[0], "Arc", e.Center.y, e.Radius, e.Center.x, "Minus"])
                        bounds.append([equationPlus, swapPoints[0], swapPoints[1], "Arc", e.Center.y, e.Radius, e.Center.x, "plus"])
                        bounds.append([equationMinus, swapPoints[1], p2, "Arc", e.Center.y, e.Radius, e.Center.x, "Minus"])
                    else:
                        bounds.append([equationPlus, p1, swapPoints[0], "Arc", e.Center.y, e.Radius, e.Center.x, "plus"])
                        bounds.append([equationMinus, swapPoints[0], swapPoints[1], "Arc", e.Center.y, e.Radius, e.Center.x, "Minus"])
                        bounds.append([equationPlus, swapPoints[1], p2, "Arc", e.Center.y, e.Radius, e.Center.x, "plus"])

                elif swap1InRange:
                    if p1.y < e.Center.y:
                        bounds.append([equationMinus, p1, swapPoints[0], "Arc", e.Center.y, e.Radius, e.Center.x, "Minus"])
                        bounds.append([equationPlus, swapPoints[0], p2, "Arc", e.Center.y, e.Radius, e.Center.x, "plus"])
                    else:
                        bounds.append([equationPlus, p1, swapPoints[0], "Arc", e.Center.y, e.Radius, e.Center.x, "plus"])
                        bounds.append([equationMinus, swapPoints[0], p2, "Arc", e.Center.y, e.Radius, e.Center.x, "Minus"])

                elif swap2InRange:
                    if p1.y < e.Center.y:
                        bounds.append([equationMinus, p1, swapPoints[1], "Arc", e.Center.y, e.Radius, e.Center.x, "Minus"])
                        bounds.append([equationPlus, swapPoints[1], p2, "Arc", e.Center.y, e.Radius, e.Center.x, "plus"])
                    else:
                        bounds.append([equationPlus, p1, swapPoints[1], "Arc", e.Center.y, e.Radius, e.Center.x, "plus"])
                        bounds.append([equationMinus, swapPoints[1], p2, "Arc", e.Center.y, e.Radius, e.Center.x, "Minus"])

                else:
                    print("ERROR")

    return bounds, edges

def areaIntegral(integrand, edges):
    bounds, edges = convertEdgesToEquations(edges)

    i = 0
    answer = 0
    errorEstimite = 0

    for b in bounds:

        if abs(b[1].x) > abs(b[2].x):
            pClosestToCenter = b[2]
            pFarthestToCenter = b[1]
        else:
            pClosestToCenter = b[1]
            pFarthestToCenter = b[2]

        mBotom = pFarthestToCenter.y / pFarthestToCenter.x
        mTop = pClosestToCenter.y / pClosestToCenter.x

        lineBotom = lambda x: mBotom * x
        lineTop = lambda x: mTop * x

        sign =  getSign(edges[i]) # needs fixxing

        if b[1].x > 0:
            if b[0] is not None:

                a, e = intg(integrand, [b[1].x, b[2].x, lineBotom, b[0], sign])
                answer += a
                errorEstimite += e

            a, e = intg(integrand, [0, b[1].x, lineBotom, lineTop, sign])
            answer += a
            errorEstimite += e
        else:
            a, e = intg(integrand, [b[2].x, 0, lineBotom, lineTop, sign])
            answer += a
            errorEstimite += e

        i += 1
    print("AREA: ", answer)

    return answer, errorEstimite

def getSign(e):
    return math.copysign(1.0, e.p2.x * e.p1.y - e.p1.x * e.p2.y) * -1

def intg(integrand, r):
    answerZone, errorEstimiteZone = integrate.dblquad(integrand, r[0], r[1], r[2], r[3])
    return abs(answerZone) * r[4], abs(errorEstimiteZone) * r[4]
