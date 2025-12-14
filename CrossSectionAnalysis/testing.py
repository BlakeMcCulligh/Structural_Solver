import math

from scipy import integrate


def triangalization(eadges, bounds):
    R = []
    i = 0
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

        sign =  getSign(eadges[i])

        if bounds[0] is not None:
            R.append([b[1].x, b[2].x, lineBotom, bounds[0], sign])

        if b[1].x > 0:
            R.append([0, b[1].x, lineBotom, lineTop, sign])
        else:
            R.append([b[2].x, 0, lineBotom, lineTop, sign])

        i += 1

def getSign(e):
    return math.copysign(1.0, e.p2.x * e.p1.y - e.p1.x * e.p2.y)

def test(y,x):
    return 1

eTop = lambda x: 0 * x + 300
eBottom = lambda x: 0.6 * x

answerZone, errorEstimiteZone = integrate.dblquad(test, 300, 500, eBottom, eTop)
print("answerZone: ", answerZone)