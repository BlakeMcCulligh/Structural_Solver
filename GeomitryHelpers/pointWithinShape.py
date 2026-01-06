import math

from Sketch.geometric_primitives.segment import Segment
from Sketch.geometric_primitives.point import Point


def point_in_shape(p, shape):

    # First: check boundary conditions
    for e in shape.Lines:
        if point_on_line_segment(p, e.p1, e.p2):
            return False
    for e in shape.Arcs:
        if point_on_arc(p, e):
            return False

    # Ray casting
    intersections = 0
    for e in shape.Lines:
        if ray_intersects_line(p, e.p1, e.p2):
            intersections += 1
    for e in shape.Arcs:
        if ray_intersects_arc(p, e):
            intersections += 1

    # Odd intersections → inside
    return intersections % 2 == 1

def ray_intersects_arc(p, arc):
    cy = arc.Center.y
    r = arc.radius
    dy = p.y - cy

    # Ray intersects circle only if |dy| <= r
    if abs(dy) > r:
        return False

    # Solve circle intersection
    dx = math.sqrt(r*r - dy*dy)

    # Intersection x-values
    xi1 = arc.Center.x - dx
    xi2 = arc.Center.x + dx

    # Only intersections to the right of point
    candidates = [x for x in (xi1, xi2) if x > p.x]

    for x_int in candidates:
        # Check if the intersection lies on the arc sweep
        px = x_int
        py = p.y
        if point_on_arc(Point(px, py), arc):
            return True

    return False

def ray_intersects_line(p, a, b):
    # Ensure a.y <= b.y
    if a.y > b.y:
        a, b = b, a

    # Does ray pass between y-values?
    if p.y < a.y or p.y >= b.y:
        return False

    # Compute x coordinate of intersection
    if almost_equal(a.y, b.y):  # horizontal edge → ignore
        return False

    x_int = a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y)

    return x_int > p.x

def almost_equal(a, b, tol=1e-9):
    return abs(a - b) < tol

def point_on_line_segment(p, a, b):
    # Check collinearity via cross product
    cross = (b.y - a.y)*(p.x - a.x) - (b.x - a.x)*(p.y - a.y)
    if abs(cross) > 1e-9:
        return False

    # Check if p is between a and b via dot product
    dot = (p.x - a.x)*(b.x - a.x) + (p.y - a.y)*(b.y - a.y)
    if dot < 0:
        return False

    squared_len = (b.x - a.x)**2 + (b.y - a.y)**2
    if dot > squared_len:
        return False

    return True


def point_on_arc(p, arc):
    # Check distance to center equals radius
    d = math.hypot(p.x - arc.Center.x, p.y - arc.Center.y)
    if not almost_equal(d, arc.radius):
        return False

    # Compute angles for arc boundaries
    ang_p1 = math.atan2(arc.p1.y - arc.Center.y, arc.p1.x - arc.Center.x)
    ang_p2 = math.atan2(arc.p2.y - arc.Center.y, arc.p2.x - arc.Center.x)
    ang_p  = math.atan2(p.y - arc.Center.y,  p.x - arc.Center.x)

    # Normalize to [0, 2π)
    def norm(a):
        while a < 0: a += 2*math.pi
        while a >= 2*math.pi: a -= 2*math.pi
        return a

    ang_p1 = norm(ang_p1)
    ang_p2 = norm(ang_p2)
    ang_p = norm(ang_p)

    # Arc always assumed CCW from N1 → N2
    if ang_p1 <= ang_p2:
        return ang_p1 - 1e-9 <= ang_p <= ang_p2 + 1e-9
    else:
        # wrap-around case
        return ang_p >= ang_p1 - 1e-9 or ang_p <= ang_p2 + 1e-9