import tkinter as tk
import numpy as np
from numpy import copy
import pyvista as pv

def triangalizeSurface(window, surface, flipNormal):

    nodes = window.PrintNodes[surface]

    cloud = pv.PolyData(nodes)
    surf = cloud.delaunay_2d()
    if flipNormal: surf.flip_faces(inplace=True)
    faces = surf.faces.reshape(-1,4)[:, 1:]
    points = surf.points

    faces = points[faces]
    return faces

def getRotationYMatrix(AngleRad):
    m = np.zeros((4,4))
    m[0,0] = np.cos(AngleRad)
    m[0,2] = np.sin(AngleRad)
    m[2,0] = -np.sin(AngleRad)
    m[1,1] = 1
    m[2,2] = np.cos(AngleRad)
    m[3,3] = 1
    return m

def getRotationXMatrix(AngleRad):
    m = np.zeros((4, 4))
    m[0,0] = 1
    m[1,1] = np.cos(AngleRad)
    m[1,2] = np.sin(AngleRad)
    m[2,1] = -np.sin(AngleRad)
    m[2,2] = np.cos(AngleRad)
    m[3,3] = 1
    return m

def getNormals(tri):
    triNormals = []
    for i in range(len(tri)):
        A = tri[i][1] - tri[i][0]
        B = tri[i][2] - tri[i][0]
        n = np.cross(A, B)
        l_norm = np.linalg.norm(n)
        n = n / l_norm
        triNormals.append(n)
    triNormals = np.array(triNormals)
    return triNormals

def removeTriFaceingAway(window, tri, triNormals):
    veiwingVector = tri[:,0] - window.Camera
    l_veiwingVector = np.linalg.norm(veiwingVector)
    veiwingVector = veiwingVector/l_veiwingVector
    result = np.array([np.dot(a, b) for a, b in zip(triNormals, veiwingVector)])
    tri = tri[result < 0]
    triNormals = triNormals[result < 0]
    return tri, triNormals

def illumination(window, solidTriNormals, numSurfTri):
    l_lightDirection = np.linalg.norm(window.LIGHT_DIR)
    window.LIGHT_DIR = window.LIGHT_DIR / l_lightDirection
    result = np.array([np.dot(a, window.LIGHT_DIR) for a in solidTriNormals])
    triColor = result * 150 + 50
    for i in range(len(triColor)):
        if triColor[i] < 50: triColor[i] = 50
        elif triColor[i] > 200: triColor[i] = 200

    triColor = triColor.tolist() + [200] * numSurfTri
    return triColor

def transformToLocal(window, node, line, tri):
    target = window.Camera + window.LookDir
    lookAtMatrix = getLookAtMatrix(window.Camera, target, window.Up)
    nodeMove = []
    for i in range(len(node)):
        nodeMove.append((np.append(node[i],1) @ lookAtMatrix)[:-1])
    nodeMove = np.array(nodeMove)
    lineMove = []
    for i in range(len(line)):
        lineSub = []
        for j in range(2):
            lineSub.append((np.append(line[i][j],1) @ lookAtMatrix)[:-1])
        lineMove.append(lineSub)
    lineMove = np.array(lineMove)
    triMove = []
    for i in range(len(tri)):
        moveSub = []
        for j in range(3):
            moveSub.append((np.append(tri[i][j], 1) @ lookAtMatrix)[:-1])
        triMove.append(moveSub)
    triMove = np.array(triMove)
    return nodeMove, lineMove, triMove

def getLookAtMatrix(postion, target, up):
    newForward = target - postion
    newForward = newForward / np.linalg.norm(newForward)
    newUp = up - (newForward * np.dot(up, newForward.T))
    newUp = newUp / np.linalg.norm(newForward)
    newRight = np.cross(newUp, newForward)
    TA = np.dot(postion,newRight)
    TB = np.dot(postion,newUp)
    TC = np.dot(postion,newForward)
    lookAtMatrix = np.array([[newRight[0],newUp[0],newForward[0],0],
                             [newRight[1],newUp[1],newForward[1],0],
                             [newRight[2],newUp[2],newForward[2],0],
                             [        -TA,     -TB,          -TC,1]])
    return lookAtMatrix

def project(window, node, line, tri):
    projectionMatrix = getProjectionMatrix(window)

    nodeProjected = []
    for i in range(len(node)):
        nodeProjected.append((np.append(node[i], 1) @ projectionMatrix) / node[i][2])
    nodeProjected = np.array(nodeProjected)

    lineProjected = []
    for i in range(len(line)):
        lineProjectedSub = []
        for j in range(2):
            lineProjectedSub.append((np.append(line[i][j],1) @ projectionMatrix) / line[i][j][2])
        lineProjected.append(lineProjectedSub)
    lineProjected = np.array(lineProjected)

    triProjected = []
    for i in range(len(tri)):
        triProjectedSub = []
        for j in range(3):
            triProjectedSub.append((np.append(tri[i][j], 1) @ projectionMatrix) / tri[i][j][2])
        triProjected.append(triProjectedSub)
    triProjected = np.array(triProjected)

    # scale to Window
    w = window.graph.winfo_width()
    h = window.graph.winfo_height()

    if len(nodeProjected) > 0:
        nodeProjected[:,[0,1]] = nodeProjected[:,[0,1]] + 1
        nodeProjected[:, 0] = nodeProjected[:, 0] * 0.5 * w
        nodeProjected[:, 1] = nodeProjected[:, 1] * 0.5 * h
    if len(lineProjected) > 0:
        lineProjected[:, :, [0, 1]] = lineProjected[:, :, [0, 1]] + 1
        lineProjected[:, :, 0] = lineProjected[:, :, 0] * 0.5 * w
        lineProjected[:, :, 1] = lineProjected[:, :, 1] * 0.5 * h
    if len(triProjected) > 0:
        triProjected[:, :, [0, 1]] = triProjected[:, :, [0, 1]] + 1
        triProjected[:, :, 0] = triProjected[:, :, 0] * 0.5 * w
        triProjected[:, :, 1] = triProjected[:, :, 1] * 0.5 * h
    return nodeProjected, lineProjected, triProjected

def getProjectionMatrix(window):
    a = window.graph.winfo_height() / window.graph.winfo_width()
    f = 1 / np.tan(np.radians(window.FOV / 2))
    q = window.Z_FAR / (window.Z_FAR - window.Z_NEAR)
    projectionMatrix = np.zeros((4, 4))
    projectionMatrix[0, 0] = a * f
    projectionMatrix[1, 1] = f
    projectionMatrix[2, 2] = q
    projectionMatrix[3, 2] = -window.Z_NEAR * q
    projectionMatrix[2, 3] = 1
    return projectionMatrix

def clipClose(node, line, tri, triColor):

    clippedNodes = []
    for i in range(len(node)):
        if node[i][2] > 0.01:
            clippedNodes.append(node[i])
    clippedNodes = np.array(clippedNodes)

    clippedLines = []
    for i in range(len(line)):
        l = Line_ClipAgainstPlane(np.array([0, 0, 0.01]), np.array([0, 0, 1]), line[i])
        if l is not None:
            clippedLines.append(l)
    clippedLines = np.array(clippedLines)

    ClippedTriangles, newColor = [], []
    for i in range(len(tri)):
        numTri, tri_out1, tri_out2 = Triangle_ClipAgainstPlane(np.array([0, 0, 0.01]), np.array([0, 0, 1]), tri[i])
        if numTri >= 1:
            ClippedTriangles.append(tri_out1)
            newColor.append(triColor[i])
        if numTri >= 2:
            ClippedTriangles.append(tri_out2)
            newColor.append(triColor[i])
    ClippedTriangles, triColor = np.array(ClippedTriangles), np.array(newColor)
    return clippedNodes, clippedLines, ClippedTriangles, triColor

def clipEadges(node, line, oldTri, oldTriColor, h, w):
    clipPlanePoint = np.array([[0, 0, 0], [0, h, 0], [0, 0, 0], [w, 0, 0]])
    clipPlaneNormal = np.array([[0, 1, 0], [0, -1, 0], [1, 0, 0], [-1, 0, 0]])

    clipedNodes = []
    for i in range(len(node)):
        if 0 < node[i][0] < w and 0 < node[i][1] < h:
            clipedNodes.append(node[i])
    clipedNodes = np.array(clipedNodes)

    for j in range(4):
        newLine = []
        for i in range(len(line)):
            l = Line_ClipAgainstPlane(clipPlanePoint[j], clipPlaneNormal[j], line[i])
            if l is not None:
                newLine.append(l)
        line = newLine
    clippedLines = np.array(line)

    for j in range(4):
        newTri, newTriColor = [], []
        for i in range(len(oldTri)):
            numTri, tri_out1, tri_out2 = Triangle_ClipAgainstPlane(clipPlanePoint[j], clipPlaneNormal[j], oldTri[i])
            if numTri >= 1:
                newTri.append(tri_out1)
                newTriColor.append(oldTriColor[i])
            if numTri >= 2:
                newTri.append(tri_out2)
                newTriColor.append(oldTriColor[i])
        oldTri, oldTriColor = newTri, newTriColor
    return clipedNodes, clippedLines, oldTri, oldTriColor

def vectorIntersectPlane(planePoint, planeNormal, lineStart, lineEnd):
    planeNormal = planeNormal/np.linalg.norm(planeNormal)
    plane_d = - np.dot(planeNormal, planePoint)
    ad, bd = np.dot(lineStart, planeNormal), np.dot(lineEnd, planeNormal)
    if bd - ad == 0: ad -= 0.0001
    t = (-plane_d - ad) / (bd - ad)
    lineStartToEnd = lineEnd - lineStart
    lineToIntersect = lineStartToEnd * t
    return lineStart + lineToIntersect

def Line_ClipAgainstPlane(planePoint, planeNormal, line):
    planeNormal = planeNormal / np.linalg.norm(planeNormal)

    def dist(p):  # return the signed shortest distance from point to plane
        p = p / np.linalg.norm(p)
        return np.dot(planeNormal, p) - np.dot(planeNormal, planePoint)

    inside_points, insidePointsCounter = [], 0
    outside_points, outsidePointsCounter = [], 0
    d = []
    for i in range(2): d.append(dist(line[i]))
    for i in range(2):
        if d[i] >= 0:
            inside_points.append(line[i])
            insidePointsCounter += 1
        else:
            outside_points.append(line[i])
            outsidePointsCounter += 1
    if insidePointsCounter == 2:
        return line
    elif insidePointsCounter == 1:
        return np.array([inside_points[0],vectorIntersectPlane(planePoint, planeNormal, inside_points[0],outside_points[0])])
    else:
        return None

def Triangle_ClipAgainstPlane(planePoint, planeNormal, in_tri):
    planeNormal = planeNormal / np.linalg.norm(planeNormal)

    def dist(p): # return the signed shortest distance from point to plane
        p = p / np.linalg.norm(p)
        return np.dot(planeNormal,p) - np.dot(planeNormal,planePoint)

    inside_points, insidePointsCounter = [], 0
    outside_points, outsidePointsCounter = [], 0
    d = []
    for i in range(3): d.append(dist(in_tri[i]))
    for i in range(3):
        if d[i] >= 0:
            inside_points.append(in_tri[i])
            insidePointsCounter += 1
        else:
            outside_points.append(in_tri[i])
            outsidePointsCounter += 1
    if insidePointsCounter == 0:
        return 0, None, None
    elif insidePointsCounter == 3:
        return 1, in_tri, None
    elif insidePointsCounter == 1:
        out_tri1 = np.array([inside_points[0],
                             vectorIntersectPlane(planePoint, planeNormal, inside_points[0], outside_points[0]),
                             vectorIntersectPlane(planePoint, planeNormal, inside_points[0], outside_points[1])])
        return 1, out_tri1, None
    elif insidePointsCounter == 2:
        out_tri1 = np.array([inside_points[0],
                             inside_points[1],
                             vectorIntersectPlane(planePoint, planeNormal, inside_points[0], outside_points[0])])
        out_tri2 = np.array([inside_points[1],
                             out_tri1[2],
                             vectorIntersectPlane(planePoint, planeNormal, inside_points[1], outside_points[0])])
        return 2, out_tri1, out_tri2
    else:
        raise Exception('ERROR')