import math
from enum import auto
from math import atan2, degrees, pi

from enum import Enum

import numpy as np
from shapely import LineString

from Sketch.geometric_primitives.common import equal_eps, v2v_angle_cw
from Sketch.geometric_primitives.point import Point, distance_p2p
from Sketch.geometric_primitives.line import Line, intersection_line_line
from Sketch.geometric_primitives.vector import Vector, cross, dot

class DIRECTION(Enum):
    CW      = auto()
    CCW     = auto()

class Arc:
    # N1, N2 -- beginning and end of the arc
    # p -- random point that belongs to the arc
    # direction is always CW (left-handed coordinate system)
    def __init__(self, p1: Point, p2: Point, p: Point):
        self.p1 = p1
        self.p2 = p2

        p1_p = Vector.from_two_points(p1, p)
        p_p2 = Vector.from_two_points(p, p2)
        p1_p2 = Vector.from_two_points(p1, p2)

        p1_p_segment_center = p1 + p1_p / 2
        p_p2_segment_center = p + p_p2 / 2
        p1_p2_segment_center = p1 + p1_p2 / 2

        l1 = Line(p1_p_segment_center, p1_p_segment_center + p1_p.rotated90ccw())
        l2 = Line(p_p2_segment_center, p_p2_segment_center + p_p2.rotated90ccw())

        center = intersection_line_line(l1, l2)

        assert center != None

        center_p1 = Vector.from_two_points(center, self.p1)
        center_p = Vector.from_two_points(center, p)

        dir = DIRECTION.CW if cross(center_p1, center_p) > 0 else DIRECTION.CCW

        if dir == DIRECTION.CCW:
            self.p1, self.p2 = self.p2, self.p1
            dir = DIRECTION.CW

        self.d = dot(Vector.from_two_points(p1_p2_segment_center, center), self.get_n())

        self.Center = None
        self.Radius = None
        self.center()
        self.radius()
        self.angle = self.angle()

    def get_n(self):
        return Vector.from_two_points(self.p1, self.p2).rotated90ccw().normalized()

    def center(self):
        p1_p2 = Vector.from_two_points(self.p1, self.p2)
        p1_p2_segment_center = self.p1 + p1_p2 / 2
        self.Center = p1_p2_segment_center + self.get_n() * self.d
        return  self.Center

    def radius(self):
        self.center()
        self.Radius = Vector.from_two_points(self.p1, self.Center).length()
        return self.Radius

    def points(self):
        return [self.p1, self.p2]

    def bb_coords(self):
        self.center()
        self.radius()
        return self.Center.x - self.Radius, self.Center.y - self.Radius, self.Center.x + self.Radius, self.Center.y + self.Radius

    def invert_direction(self):
        self.p1, self.p2 = self.p2, self.p1

    def middle_point(self):
        self.center()

        c_p1 = Vector.from_two_points(self.Center, self.p1)
        c_p2 = Vector.from_two_points(self.Center, self.p2)

        angle = v2v_angle_cw(c_p1, c_p2)

        return self.Center + c_p1.rotated(angle / 2)

    def angle(self):
        self.center()
        v1 = Vector.from_two_points(self.Center, self.p1)
        v2 = Vector.from_two_points(self.Center, self.p2)

        return v2v_angle_cw(v1, v2)

def distance_p2a(p: Point, arc: Arc):
    arc_center = arc.center()

    c_p1 = Vector.from_two_points(arc_center, arc.p1)
    c_p2 = Vector.from_two_points(arc_center, arc.p2)
    c_p = Vector.from_two_points(arc_center, p)

    is_inside_sector = equal_eps(v2v_angle_cw(c_p1, c_p) + v2v_angle_cw(c_p, c_p2), v2v_angle_cw(c_p1, c_p2))

    if is_inside_sector:
        return abs(arc.radius() - distance_p2p(arc_center, p))
    else:
        return min(distance_p2p(p, arc.p1), distance_p2p(p, arc.p2))

def create_arc_segment(center_x, center_y, radius, start_angle_deg, end_angle_deg, num_segments=100):
    """
    Generates a shapely LineString approximating an arc segment.

    :param center_x: X coordinate of the center.
    :param center_y: Y coordinate of the center.
    :param radius: The radius of the arc.
    :param start_angle_deg: The start angle in degrees (e.g., 0).
    :param end_angle_deg: The end angle in degrees (e.g., 90).
    :param num_segments: The number of linear segments to use for approximation.
    :return: A shapely LineString geometry.
    """
    # Convert angles from degrees to radians
    start_angle_rad = math.radians(start_angle_deg)
    end_angle_rad = math.radians(end_angle_deg)

    # Generate angles for the segments
    # The order of start/end angles may need adjustment depending on desired arc direction (CW/CCW)
    angles = np.linspace(start_angle_rad, end_angle_rad, num_segments)

    # Calculate coordinates
    x = center_x + radius * np.cos(angles)
    y = center_y + radius * np.sin(angles)

    # Create the LineString from the coordinates
    arc_coords = np.column_stack([x, y])
    arc_segment = LineString(arc_coords)

    return arc_segment


def get_angle_of_point(center_x, center_y, point_x, point_y):
    """
    Calculates the angle of a point on a circle relative to the center.

    The angle is measured counter-clockwise from the positive x-axis (3 o'clock position).
    Returns the angle in degrees in the range [-180, 180].
    """
    # Calculate the difference in coordinates
    dx = point_x - center_x
    dy = point_y - center_y

    # Use math.atan2 to get the angle in radians
    angle_radians = math.atan2(dy, dx)

    # Convert radians to degrees
    angle_degrees = math.degrees(angle_radians)

    return angle_degrees


    
