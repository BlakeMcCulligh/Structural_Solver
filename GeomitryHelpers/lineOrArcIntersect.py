import math

from Sketch.geometric_primitives.arc import Arc
from Sketch.geometric_primitives.point import Point
from Sketch.geometric_primitives.segment import Segment



def checkIfLineOrArcIntersect(e1, e2):
    if isinstance(e1, Segment) and isinstance(e2, Segment):
        return line_line_intersect(e1, e2)

    if isinstance(e1, Segment) and isinstance(e2, Arc):
        return line_arc_intersect(e1, e2)

    if isinstance(e1, Arc) and isinstance(e2, Segment):
        return line_arc_intersect(e2, e1)

    if isinstance(e1, Arc) and isinstance(e2, Arc):
        return arc_arc_intersect(e1, e2)

    return False

def line_line_intersect(l1, l2):
    p = l1.p1
    r = Point(l1.p2.x - l1.p1.x, l1.p2.y - l1.p1.y)

    q = l2.p1
    s = Point(l2.p2.x - l2.p1.x, l2.p2.y - l2.p1.y)

    def cross(a, b):
        return a.x * b.y - a.y * b.x

    rxs = cross(r, s)
    q_p = Point(q.x - p.x, q.y - p.y)
    q_pxr = cross(q_p, r)

    # Collinear or parallel ⇒ no interior intersection
    if abs(rxs) < 1e-12:
        return False

    t = cross(q_p, s) / rxs
    u = cross(q_p, r) / rxs

    # Strict interior intersection
    if 1e-12 < t < 1 - 1e-12 and 1e-12 < u < 1 - 1e-12:
        return True

    return False

def point_on_arc_strict(p, arc):
    # Distance from Center must equal radius
    d = math.hypot(p.x - arc.Center.x, p.y - arc.Center.y)
    if not almost_equal(d, arc.radius):
        return False

    # Compute angles
    def norm(a):
        while a < 0:
            a += 2*math.pi
        while a >= 2*math.pi:
            a -= 2*math.pi
        return a

    ang1 = norm(math.atan2(arc.p1.y - arc.Center.y, arc.p1.x - arc.Center.x))
    ang2 = norm(math.atan2(arc.p2.y - arc.Center.y, arc.p2.x - arc.Center.x))
    angp = norm(math.atan2(p.y - arc.Center.y, p.x - arc.Center.x))

    # Check interior of arc (CCW assumed)
    if ang1 <= ang2:
        return ang1 < angp < ang2
    else:
        return angp > ang1 or angp < ang2


def line_arc_intersect(line, arc):
    # Shift line so arc.Center becomes (0,0)
    x1 = line.p1.x - arc.Center.x
    y1 = line.p1.y - arc.Center.y
    x2 = line.p2.x - arc.Center.x
    y2 = line.p2.y - arc.Center.y

    dx = x2 - x1
    dy = y2 - y1

    # Quadratic coefficients for intersection with circle
    a = dx * dx + dy * dy
    b = 2 * (x1 * dx + y1 * dy)
    c = x1 * x1 + y1 * y1 - arc.radius * arc.radius

    disc = b * b - 4 * a * c
    if disc < 0:
        return False

    sqrt_disc = math.sqrt(max(disc, 0))

    # Up to two intersection values
    t_candidates = [
        (-b - sqrt_disc) / (2 * a),
        (-b + sqrt_disc) / (2 * a)
    ]

    for t in t_candidates:
        # Strictly within segment
        if t <= 1e-12 or t >= 1 - 1e-12:
            continue

        px = line.p1.x + t * (line.p2.x - line.p1.x)
        py = line.p1.y + t * (line.p2.y - line.p1.y)
        ip = Point(px, py)

        if point_on_arc_strict(ip, arc):
            return True

    return False

def arc_arc_intersect(a1, a2):
    cx1, cy1 = a1.Center.x, a1.Center.y
    cx2, cy2 = a2.Center.x, a2.Center.y
    r1, r2 = a1.radius, a2.radius

    dx = cx2 - cx1
    dy = cy2 - cy1
    d = math.hypot(dx, dy)

    # No intersection
    if d > r1 + r2 or d < abs(r1 - r2) or d == 0:
        return False

    # Distance from Center1 to intersection line
    a = (r1*r1 - r2*r2 + d*d) / (2*d)
    h = math.sqrt(max(r1*r1 - a*a, 0))

    xm = cx1 + a * dx / d
    ym = cy1 + a * dy / d

    rx = -dy * (h/d)
    ry = dx * (h/d)

    candidates = [
        Point(xm + rx, ym + ry),
        Point(xm - rx, ym - ry)
    ]

    for p in candidates:
        if point_on_arc_strict(p, a1) and point_on_arc_strict(p, a2):
            return True

    return False

def almost_equal(a, b, tol=1e-9):
    return abs(a - b) < tol