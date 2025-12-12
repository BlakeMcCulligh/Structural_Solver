# helpers

from collections import Counter
from Sketch.geometric_primitives.arc import Arc
from Sketch.geometric_primitives.line import Line, distance_p2l
from Sketch.geometric_primitives.point import distance_p2p, Point
from Sketch.geometric_primitives.segment import Segment
from Sketch.geometric_primitives.vector import Vector, cross, dot

def pairs(entities):
    return zip(entities[:-1], entities[1:])

def all_entities(entities, cls):
    return all(entity.__class__ is cls for entity in entities)

# constraints

def parallel(s1: Segment, s2: Segment):
    return [cross(s1.vector(), s2.vector())]

def perpendicular(s1: Segment, s2: Segment):
    return [dot(s1.vector(), s2.vector())]

def equal_length_or_radius(entity1, entity2):
    entities = (entity1, entity2)
    segments, arcs = all_entities(entities, Segment), all_entities(entities, Arc)
    assert segments or arcs
    return [entity1.length() - entity2.length()] if segments else [entity1.radius - entity2.radius]

def length_or_radius_equal_to_Number(entity1, entity2, entity3 = None):

    if entity3 is None:
        if isinstance(entity1, float):
            return length(entity2, entity1)
        else:
            return length(entity1, entity2)
    else:
        if isinstance(entity3, float):
            return [(Vector(entity2.x - entity1.x, entity2.y - entity1.y).length() - entity3)]
        else:
            return [(Vector(entity3.x - entity2.x, entity3.y - entity2.y).length() - entity1)]

def length(segment: Segment, length: float):
    return [(Vector(segment.p2.x - segment.p1.x, segment.p2.y - segment.p1.y).length() - length)]

def tangency(entity1, entity2):
    temp = Counter((entity1.__class__, entity2.__class__))
    if temp == Counter((Arc, Segment)):
        arc = entity1 if isinstance(entity1, Arc) else entity2
        segment = entity1 if isinstance(entity1, Segment) else entity2
        return [distance_p2l(arc.center(), Line(segment.p1, segment.p2)) - arc.radius]
    elif temp == Counter((Arc, Arc)):
        return [distance_p2p(entity1.center(), entity2.center()) - (entity1.radius() + entity2.radius())]
    else:
        assert False

def concentricity(arc1: Arc, arc2: Arc):
    return [distance_p2p(arc1.center(), arc2.center())]