"""
Handels the calculations of converting the 3D space to be printed on the screen. Credit to onelonecoder/javidx9 for
tutorial that can be found at https://youtu.be/ih20l3pJoeU.
"""

import numpy as np
import pyvista as pv

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh","onelonecoder/javidx9"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

def triangalize_surface(window, surface, flip_normal):
    """
    Splits a sirface into multiple triangular surfaces.

    :param window: Object storing the main window.
    :param surface: List of nodes that form the surface to be split up.
    :param flip_normal: If the normals of the new surface should be in the negative direction.
    :return: List of nodes, that for triangular surfaces.  Shape: (# Triangular Surfaces, 3, 3)
    """

    nodes = window.PrintNodes[surface]

    cloud = pv.PolyData(nodes)
    surf = cloud.delaunay_2d()
    if flip_normal: surf.flip_faces(inplace=True)
    faces = surf.faces.reshape(-1,4)[:, 1:]
    points = surf.points

    faces = points[faces]
    return faces

def get_rotation_y_matrix(angle_rad):
    """
    Gets the Y rotation matrix for the given angle in radians.

    :param angle_rad: Angle in radians to get the rotaton matrix for.
    :return: Y rotation matrix.
    """

    m = np.zeros((4,4))
    m[0,0] = np.cos(angle_rad)
    m[0,2] = np.sin(angle_rad)
    m[2,0] = -np.sin(angle_rad)
    m[1,1] = 1
    m[2,2] = np.cos(angle_rad)
    m[3,3] = 1
    return m

def get_rotation_x_matrix(angle_rad):
    """
    Gets the X rotation matrix for the given angle in radians.

    :param angle_rad: Angle in radians to get the rotaton matrix for.
    :return: X rotation matrix.
    """

    m = np.zeros((4, 4))
    m[0,0] = 1
    m[1,1] = np.cos(angle_rad)
    m[1,2] = np.sin(angle_rad)
    m[2,1] = -np.sin(angle_rad)
    m[2,2] = np.cos(angle_rad)
    m[3,3] = 1
    return m

def get_normals(tri):
    """
    Gets the normals for triangular surfaces.

    :param tri: List of nodes that forms the triangular surfaces to get the normals for.
                Shape: (# Triangular Surfaces, 3, 3)
    :return: List of vectors normal to the surfaces.
    """

    tri_normals = []
    for i in range(len(tri)):
        a = tri[i][1] - tri[i][0]
        b = tri[i][2] - tri[i][0]
        n = np.cross(a, b)
        l_norm = np.linalg.norm(n)
        n = n / l_norm
        tri_normals.append(n)
    tri_normals = np.array(tri_normals)
    return tri_normals

def remove_tri_faceing_away(window, tri, tri_normals):
    """
    Removes all triangular surfaces that are facing away from the camera.

    :param window: Object storing the main window.
    :param tri: List of nodes that forms the triangular surfaces. Shape: (# Triangular Surfaces, 3, 3)
    :param tri_normals: List of vectors normal to the triangular surfaces.
    :return: List of nodes that forms the triangular surfaces, List of vectors normal to the triangular surfaces.
    """

    veiwing_vector = tri[:,0] - window.Camera
    l_veiwing_vector = np.linalg.norm(veiwing_vector)
    veiwing_vector = veiwing_vector/l_veiwing_vector
    result = np.array([np.dot(a, b) for a, b in zip(tri_normals, veiwing_vector)])
    tri = tri[result < 0]
    tri_normals = tri_normals[result < 0]
    return tri, tri_normals

def illumination(window, solid_tri_normals, num_surf_tri):
    """
    Gets the color variable for each triangular surface. The color variable is an int
    between 0 and 255 and is used for all 3 values in rgb.

    :param window: Object storing the main window.
    :param solid_tri_normals: List of vectors normal to the solid's triangular surfaces.
    :param num_surf_tri: Number of surface triangular surfaces.
    :return: Color Variable between 0 and 255.
    """

    l_light_direction = np.linalg.norm(window.LIGHT_DIR)
    window.LIGHT_DIR = window.LIGHT_DIR / l_light_direction
    result = np.array([np.dot(a, window.LIGHT_DIR) for a in solid_tri_normals])
    tri_color = result * 150 + 50
    for i in range(len(tri_color)):
        if tri_color[i] < 50: tri_color[i] = 50
        elif tri_color[i] > 200: tri_color[i] = 200

    tri_color = tri_color.tolist() + [200] * num_surf_tri
    return tri_color

def transform_to_local(window, node, line, tri):
    """
    Transforms the nodes, lines and surfaces to local coordinates relative to the camera.

    :param window: Object storing the main window.
    :param node: List of nodes. Shape: (# Nodes, 3)
    :param line: List of nodes that forms lines. Shape: (# Lines, 2, 3)
    :param tri: List of nodes that forms triagular surfaces. Shape: (# Triangular Surfaces, 3, 3)
    :return: Nodes, Lines, Surfaces.
    """

    target = window.Camera + window.LookDir
    look_at_matrix = get_look_at_matrix(window.Camera, target, window.Up)
    node_move = []
    for i in range(len(node)):
        node_move.append((np.append(node[i],1) @ look_at_matrix)[:-1])
    node_move = np.array(node_move)
    line_move = []
    for i in range(len(line)):
        line_sub = []
        for j in range(2):
            line_sub.append((np.append(line[i][j],1) @ look_at_matrix)[:-1])
        line_move.append(line_sub)
    line_move = np.array(line_move)
    tri_move = []
    for i in range(len(tri)):
        move_sub = []
        for j in range(3):
            move_sub.append((np.append(tri[i][j], 1) @ look_at_matrix)[:-1])
        tri_move.append(move_sub)
    tri_move = np.array(tri_move)
    return node_move, line_move, tri_move

def get_look_at_matrix(postion, target, up):
    """
    Gets the look at matrix for the given position, toarget, and up direction of the camera.

    :param postion: Position vector of the camera.
    :param target: Position vector of the coameras target.
    :param up: Vector pointing in the up direction of the camera.
    :return: Look at matrix for the camera.
    """

    new_forward = target - postion
    new_forward = new_forward / np.linalg.norm(new_forward)
    new_up = up - (new_forward * np.dot(up, new_forward.T))
    new_up = new_up / np.linalg.norm(new_forward)
    new_right = np.cross(new_up, new_forward)
    TA = np.dot(postion,new_right)
    TB = np.dot(postion,new_up)
    TC = np.dot(postion,new_forward)
    look_at_matrix = np.array([[new_right[0],new_up[0],new_forward[0],0],
                             [new_right[1],new_up[1],new_forward[1],0],
                             [new_right[2],new_up[2],new_forward[2],0],
                             [        -TA,     -TB,          -TC,1]])
    return look_at_matrix

def project(window, node, line, tri):
    """
    Projects the nodes, lines, and surfaces and scales them to the size of the window.

    :param window: Object storing the main window.
    :param node: List of nodes. Shape: (# Nodes, 3)
    :param line: List of nodes that forms lines. Shape: (# Lines, 2, 3)
    :param tri: List of nodes that forms triagular surfaces. Shape: (# Triangular Surfaces, 3, 3)
    :return: Nodes, Lines, Surfaces.
    """

    projection_matrix = get_projection_matrix(window)

    node_projected = []
    for i in range(len(node)):
        node_projected.append((np.append(node[i], 1) @ projection_matrix) / node[i][2])
    node_projected = np.array(node_projected)

    line_projected = []
    for i in range(len(line)):
        line_projected_sub = []
        for j in range(2):
            line_projected_sub.append((np.append(line[i][j],1) @ projection_matrix) / line[i][j][2])
        line_projected.append(line_projected_sub)
    line_projected = np.array(line_projected)

    tri_projected = []
    for i in range(len(tri)):
        tri_projected_sub = []
        for j in range(3):
            tri_projected_sub.append((np.append(tri[i][j], 1) @ projection_matrix) / tri[i][j][2])
        tri_projected.append(tri_projected_sub)
    tri_projected = np.array(tri_projected)

    # scale to Window
    w = window.graph.winfo_width()
    h = window.graph.winfo_height()

    if len(node_projected) > 0:
        node_projected[:,[0,1]] = node_projected[:,[0,1]] + 1
        node_projected[:, 0] = node_projected[:, 0] * 0.5 * w
        node_projected[:, 1] = node_projected[:, 1] * 0.5 * h
    if len(line_projected) > 0:
        line_projected[:, :, [0, 1]] = line_projected[:, :, [0, 1]] + 1
        line_projected[:, :, 0] = line_projected[:, :, 0] * 0.5 * w
        line_projected[:, :, 1] = line_projected[:, :, 1] * 0.5 * h
    if len(tri_projected) > 0:
        tri_projected[:, :, [0, 1]] = tri_projected[:, :, [0, 1]] + 1
        tri_projected[:, :, 0] = tri_projected[:, :, 0] * 0.5 * w
        tri_projected[:, :, 1] = tri_projected[:, :, 1] * 0.5 * h
    return node_projected, line_projected, tri_projected

def get_projection_matrix(window):
    """
    Constructs the projection matrix for the given camera properties.

    :param window: Object storing the main window.
    :return: Projection Matrix
    """

    a = window.graph.winfo_height() / window.graph.winfo_width()
    f = 1 / np.tan(np.radians(window.FOV / 2))
    q = window.Z_FAR / (window.Z_FAR - window.Z_NEAR)
    projection_matrix = np.zeros((4, 4))
    projection_matrix[0, 0] = a * f
    projection_matrix[1, 1] = f
    projection_matrix[2, 2] = q
    projection_matrix[3, 2] = -window.Z_NEAR * q
    projection_matrix[2, 3] = 1
    return projection_matrix

def clip_close(node, line, tri, tri_color):
    """
    Cuts nodes, lines, and surfaces that are to close or behined the camera.

    :param node: List of nodes. Shape: (# Nodes, 3)
    :param line: List of nodes that forms lines. Shape: (# Lines, 2, 3)
    :param tri: List of nodes that forms triagular surfaces. Shape: (# Triangular Surfaces, 3, 3)
    :param tri_color: Color Variable between 0 and 255.
    :return: Nodes, Lines, Triangles, Triangle Colors.
    """

    clipped_nodes = []
    for i in range(len(node)):
        if node[i][2] > 0.01:
            clipped_nodes.append(node[i])
    clipped_nodes = np.array(clipped_nodes)

    clipped_lines = []
    for i in range(len(line)):
        l = line_clip_against_plane(np.array([0, 0, 0.01]), np.array([0, 0, 1]), line[i])
        if l is not None:
            clipped_lines.append(l)
    clipped_lines = np.array(clipped_lines)

    clipped_triangles, new_color = [], []
    for i in range(len(tri)):
        num_tri, tri_out1, tri_out2 = triangle_clip_against_plane(np.array([0, 0, 0.01]), np.array([0, 0, 1]), tri[i])
        if num_tri >= 1:
            clipped_triangles.append(tri_out1)
            new_color.append(tri_color[i])
        if num_tri >= 2:
            clipped_triangles.append(tri_out2)
            new_color.append(tri_color[i])
    clipped_triangles, tri_color = np.array(clipped_triangles), np.array(new_color)
    return clipped_nodes, clipped_lines, clipped_triangles, tri_color

def clip_eadges(node, line, old_tri, old_tri_color, h, w):
    """
    Cuts nodes, lines, and surfaces that are ouside the eadges of the camera.

    :param node:  List of nodes. Shape: (# Nodes, 3)
    :param line:  List of nodes that forms lines. Shape: (# Lines, 2, 3)
    :param old_tri:  List of nodes that forms triagular surfaces. Shape: (# Triangular Surfaces, 3, 3)
    :param old_tri_color:  Color Variable between 0 and 255.
    :param h: Height of window.
    :param w: Width of window.
    :return: Nodes, Lines, Triangles, Triangle Colors.
    """

    clip_plane_point = np.array([[0, 0, 0], [0, h, 0], [0, 0, 0], [w, 0, 0]])
    clip_plane_normal = np.array([[0, 1, 0], [0, -1, 0], [1, 0, 0], [-1, 0, 0]])

    cliped_nodes = []
    for i in range(len(node)):
        if 0 < node[i][0] < w and 0 < node[i][1] < h:
            cliped_nodes.append(node[i])
    cliped_nodes = np.array(cliped_nodes)

    for j in range(4):
        new_line = []
        for i in range(len(line)):
            l = line_clip_against_plane(clip_plane_point[j], clip_plane_normal[j], line[i])
            if l is not None:
                new_line.append(l)
        line = new_line
    clipped_lines = np.array(line)

    for j in range(4):
        new_tri, new_tri_color = [], []
        for i in range(len(old_tri)):
            num_tri, tri_out1, tri_out2 = triangle_clip_against_plane(clip_plane_point[j], clip_plane_normal[j], old_tri[i])
            if num_tri >= 1:
                new_tri.append(tri_out1)
                new_tri_color.append(old_tri_color[i])
            if num_tri >= 2:
                new_tri.append(tri_out2)
                new_tri_color.append(old_tri_color[i])
        old_tri, old_tri_color = new_tri, new_tri_color
    return cliped_nodes, clipped_lines, old_tri, old_tri_color

def vector_intersect_plane(plane_point, plane_normal, line_start, line_end):
    """
    Fineds the location of where a line and a plane intersect.

    :param plane_point: Node on the plane.
    :param plane_normal: Normal vector of the plane.
    :param line_start: Start node of the line.
    :param line_end: End node of the line.
    :return: Node where the plane and line intersect.
    """

    plane_normal = plane_normal / np.linalg.norm(plane_normal)
    plane_d = - np.dot(plane_normal, plane_point)
    ad, bd = np.dot(line_start, plane_normal), np.dot(line_end, plane_normal)
    if bd - ad == 0: ad -= 0.0001
    t = (-plane_d - ad) / (bd - ad)
    line_start_to_end = line_end - line_start
    line_to_intersect = line_start_to_end * t
    return line_start + line_to_intersect

def line_clip_against_plane(plane_point, plane_normal, line):
    """
    Cuts a line where it intersects with a plane.

    :param plane_point: Node onn the plane.
    :param plane_normal: Normal vector of the plane.
    :param line: Start and end nodes of the line.
    :return: Start and end nodes of the line on the inside of the plane.
    """

    plane_normal = plane_normal / np.linalg.norm(plane_normal)

    def dist(p):  # return the signed shortest distance from point to plane
        p = p / np.linalg.norm(p)
        return np.dot(plane_normal, p) - np.dot(plane_normal, plane_point)

    inside_points, inside_points_counter = [], 0
    outside_points, outside_points_counter = [], 0
    d = []
    for i in range(2): d.append(dist(line[i]))
    for i in range(2):
        if d[i] >= 0:
            inside_points.append(line[i])
            inside_points_counter += 1
        else:
            outside_points.append(line[i])
            outside_points_counter += 1
    if inside_points_counter == 2:
        return line
    elif inside_points_counter == 1:
        return np.array([inside_points[0], vector_intersect_plane(plane_point, plane_normal, inside_points[0], outside_points[0])])
    else:
        return None

def triangle_clip_against_plane(plane_point, plane_normal, in_tri):
    """
    Cuts a triangle where it intersects with a plane.

    :param plane_point: Point on the plane.
    :param plane_normal: Normal vector of the plane.
    :param in_tri: Three Corner nodes of the triangle.
    :return: Num Triangles, Tri 1, Tri 2 (Tris can be none depending on the number of triangles).
    """

    plane_normal = plane_normal / np.linalg.norm(plane_normal)

    def dist(p): # return the signed shortest distance from point to plane
        p = p / np.linalg.norm(p)
        return np.dot(plane_normal, p) - np.dot(plane_normal, plane_point)

    inside_points, inside_points_counter = [], 0
    outside_points, outside_points_counter = [], 0
    d = []
    for i in range(3): d.append(dist(in_tri[i]))
    for i in range(3):
        if d[i] >= 0:
            inside_points.append(in_tri[i])
            inside_points_counter += 1
        else:
            outside_points.append(in_tri[i])
            outside_points_counter += 1
    if inside_points_counter == 0:
        return 0, None, None
    elif inside_points_counter == 3:
        return 1, in_tri, None
    elif inside_points_counter == 1:
        out_tri1 = np.array([inside_points[0],
                             vector_intersect_plane(plane_point, plane_normal, inside_points[0], outside_points[0]),
                             vector_intersect_plane(plane_point, plane_normal, inside_points[0], outside_points[1])])
        return 1, out_tri1, None
    elif inside_points_counter == 2:
        out_tri1 = np.array([inside_points[0],
                             inside_points[1],
                             vector_intersect_plane(plane_point, plane_normal, inside_points[0], outside_points[0])])
        out_tri2 = np.array([inside_points[1],
                             out_tri1[2],
                             vector_intersect_plane(plane_point, plane_normal, inside_points[1], outside_points[0])])
        return 2, out_tri1, out_tri2
    else:
        raise Exception('ERROR')