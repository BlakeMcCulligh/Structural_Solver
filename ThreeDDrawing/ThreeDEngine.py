import time
import tkinter as tk
import keyboard
import numpy as np
from numpy import copy

class MainWindow(tk.Frame):
    def __init__(self, root, triangles):
        super().__init__()
        self.sizeW, self.sizeH = 1000, 800
        self.root = root
        self.root.geometry(str(self.sizeW) + "x" + str(self.sizeH))
        self.graph = tk.Canvas(root, bg="black")
        self.graph.pack(fill="both", expand=True)
        self.camera = np.array([0.0,0.0,0.0])
        self.up = np.array([0.0,1.0,0.0])
        self.lookDir = np.array([0.0,0.0,1.0])
        self.cam_Y_rot = 0
        self.cam_X_rot = 0
        self.target = np.array([0.0,0.0,1.0])
        self.FOV = 90
        self.Z_far, self.Z_near = 1000, 0.1
        self.light_direction = np.array([0,0,-1])

        self.tri = triangles

        self.root.bind("<MouseWheel>", self.zoom)

        self.scroolCordsLast = []
        self.shiftStateLast = False
        self.root.bind("<Button-2>", self.scroolDown)
        self.root.bind("<B2-Motion>", self.scroolUpdate)
        self.root.bind("<ButtonRelease-2>", self.scroolUp)

    def scroolDown(self, event):
        if bool(event.state & 0x0001):
            self.shiftStateLast = True
            self.rotStart(event)
        else:
            self.shiftStateLast = False
            self.panStart(event)

    def scroolUpdate(self, event):
        if bool(event.state & 0x0001):
            if not self.shiftStateLast:
                self.shiftStateLast = True
                self.panEnd(event)
                self.rotStart(event)
            else: self.rotUpdate(event)
        else:
            if self.shiftStateLast:
                self.shiftStateLast = False
                self.rotEnd(event)
                self.panStart(event)
            else: self.panUpdate(event)

    def scroolUp(self, event):
        if bool(event.state & 0x0001):
            if not self.shiftStateLast: self.panEnd(event)
            else: self.rotEnd(event)
        else:
            if self.shiftStateLast: self.rotEnd(event)
            else: self.panEnd(event)

    def zoom(self, event):
        if event.delta > 0:
            self.camera = self.camera + self.lookDir * 0.1
        else:
            self.camera = self.camera - self.lookDir * 0.1
        self.updateTriangles()

    def panStart(self, event):
        self.scroolCordsLast = [event.x, event.y]

    def panUpdate(self, event):
        cords = [event.x, event.y]
        self.panMove(cords)

    def panEnd(self, event):
        cords = [event.x, event.y]
        self.panMove(cords)

    def panMove(self, cords):
        changeCords = [cords[0] - self.scroolCordsLast[0], cords[1] - self.scroolCordsLast[1]]
        self.camera = self.camera - self.up * changeCords[1] * 0.01
        self.camera = self.camera + np.cross(self.lookDir, self.up) * changeCords[0] * 0.01
        self.updateTriangles()
        self.scroolCordsLast = cords

    def rotStart(self, event):
        self.scroolCordsLast = [event.x, event.y]

    def rotUpdate(self, event):
        cords = [event.x, event.y]
        self.rotMove(cords)

    def rotEnd(self, event):
        cords = [event.x, event.y]
        self.rotMove(cords)

    def rotMove(self, cords):
        changeCords = [cords[0] - self.scroolCordsLast[0], cords[1] - self.scroolCordsLast[1]]
        self.cam_Y_rot = self.cam_Y_rot + changeCords[0] * 0.005
        self.cam_X_rot = self.cam_X_rot + changeCords[1] * 0.005
        self.updateTriangles()
        self.scroolCordsLast = cords

    def updateTriangles(self):
        tri = copy(self.tri)
        self.updateCamData()
        triNormals = getNormals(tri)
        tri, triNormals = removeTriFaceingAway(self, tri, triNormals)
        triColor = illumination(self, triNormals)
        tri = transformToLocal(self, tri)
        tri, triColor = clipTriClose(tri, triColor)
        self.graph.delete("all")
        if len(tri) > 0:
            w = self.graph.winfo_width()
            h = self.graph.winfo_height()
            tri = projectTri(self, tri)
            tri, triColor = sortPrintingOrder(tri, triColor)
            oldTri, oldTriColor = clipTriEadges(tri[:,:,:-1], triColor, h, w)
            # print triangles
            triPrint = []
            for i in range(len(oldTri)):
                triPrint.append([oldTri[i][0][0],oldTri[i][0][1],
                                     oldTri[i][1][0],oldTri[i][1][1],
                                     oldTri[i][2][0],oldTri[i][2][1]])
            for i in range(len(triPrint)):
                color = int(oldTriColor[i])
                self.graph.create_polygon(triPrint[i],outline='',fill="#%02x%02x%02x" % (color, color, color))
        self.update()

    def updateCamData(self):
        matCameraRotY = getRotationYMatrix(self.cam_Y_rot)
        matCameraRotX = getRotationXMatrix(self.cam_X_rot)
        self.lookDir = (np.append(np.array([0, 0, 1]), 1) @ matCameraRotY @ matCameraRotX)[:-1]
        self.lookDir = self.lookDir / np.linalg.norm(self.lookDir)
        self.target = self.camera + self.lookDir
        self.up = (np.append(np.array([0, 1, 0]), 1) @ matCameraRotY @ matCameraRotX)[:-1]
        self.up = self.up / np.linalg.norm(self.up)

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
    veiwingVector = tri[:,0] - window.camera
    l_veiwingVector = np.linalg.norm(veiwingVector)
    veiwingVector = veiwingVector/l_veiwingVector
    result = np.array([np.dot(a, b) for a, b in zip(triNormals, veiwingVector)])
    tri = tri[result < 0]
    triNormals = triNormals[result < 0]
    return tri, triNormals

def illumination(window, triNormals):
    l_lightDirection = np.linalg.norm(window.light_direction)
    window.light_direction = window.light_direction / l_lightDirection
    result = np.array([np.dot(a, window.light_direction) for a in triNormals])
    triColor = result * 255
    for i in range(len(triColor)):
        if triColor[i] < 50: triColor[i] = 50
        elif triColor[i] > 255: triColor[i] = 255
    return triColor

def transformToLocal(window, tri):
    target = window.camera + window.lookDir
    lookAtMatrix = getLookAtMatrix(window.camera, target, window.up)
    triMove = []
    for i in range(len(tri)):
        moveSub = []
        for j in range(3):
            moveSub.append((np.append(tri[i][j], 1) @ lookAtMatrix)[:-1])
        triMove.append(moveSub)
    triMove = np.array(triMove)
    return triMove

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

def projectTri(window, tri):
    projectionMatrix = getProjectionMatrix(window)
    triProjected = []
    for i in range(len(tri)):
        triProjectedSub = []
        for j in range(3):
            triProjectedSub.append((np.append(tri[i][j], 1) @ projectionMatrix) / tri[i][j][2])
        triProjected.append(triProjectedSub)
    triProjected = np.array(triProjected)
    # scale to window
    triProjected[:, :, [0, 1]] = triProjected[:, :, [0, 1]] + 1
    triProjected[:, :, 0] = triProjected[:, :, 0] * 0.5 * window.graph.winfo_width()
    triProjected[:, :, 1] = triProjected[:, :, 1] * 0.5 * window.graph.winfo_height()
    return triProjected

def getProjectionMatrix(window):
    a = window.graph.winfo_height() / window.graph.winfo_width()
    f = 1 / np.tan(np.radians(window.FOV / 2))
    q = window.Z_far / (window.Z_far - window.Z_near)
    projectionMatrix = np.zeros((4, 4))
    projectionMatrix[0, 0] = a * f
    projectionMatrix[1, 1] = f
    projectionMatrix[2, 2] = q
    projectionMatrix[3, 2] = -window.Z_near * q
    projectionMatrix[2, 3] = 1
    return projectionMatrix

def clipTriClose(tri, triColor):
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
    return ClippedTriangles, triColor

def clipTriEadges(oldTri, oldTriColor, h, w):
    clipPlanePoint = np.array([[0, 0, 0], [0, h, 0], [0, 0, 0], [w, 0, 0]])
    clipPlaneNormal = np.array([[0, 1, 0], [0, -1, 0], [1, 0, 0], [-1, 0, 0]])
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
    return oldTri, oldTriColor

def vectorIntersectPlane(planePoint, planeNormal, lineStart, lineEnd):
    planeNormal = planeNormal/np.linalg.norm(planeNormal)
    plane_d = - np.dot(planeNormal, planePoint)
    ad, bd = np.dot(lineStart, planeNormal), np.dot(lineEnd, planeNormal)
    if bd - ad == 0: ad -= 0.0001
    t = (-plane_d - ad) / (bd - ad)
    lineStartToEnd = lineEnd - lineStart
    lineToIntersect = lineStartToEnd * t
    return lineStart + lineToIntersect

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

def sortPrintingOrder(tri, triColor):
    avgDistance = np.sum(tri[:, :, 2], axis=-1)
    idx = np.argsort(avgDistance)
    tri = tri[idx[::-1]]
    triColor = triColor[idx[::-1]]
    return tri, triColor

def main():

    # Window Set Up
    root_widget = tk.Tk()
    root_widget.title('3D Grphics Engine')

    # Geomitry Set Up
    nodes: np.ndarray = np.empty((0, 3))
    triangles: np.ndarray = np.empty((0, 3, 3))

    # Defining Cube
    nodes = np.vstack((nodes, [0, 0, 0]))
    nodes = np.vstack((nodes, [1, 0, 0]))
    nodes = np.vstack((nodes, [1, 1, 0]))
    nodes = np.vstack((nodes, [1, 0, 1]))
    nodes = np.vstack((nodes, [1, 1, 1]))
    nodes = np.vstack((nodes, [0, 1, 0]))
    nodes = np.vstack((nodes, [0, 1, 1]))
    nodes = np.vstack((nodes, [0, 0, 1]))
    triangles = np.vstack((triangles, [[nodes[0],nodes[5],nodes[2]]]))
    triangles = np.vstack((triangles, [[nodes[0],nodes[2],nodes[1]]]))
    triangles = np.vstack((triangles, [[nodes[1], nodes[2], nodes[4]]]))
    triangles = np.vstack((triangles, [[nodes[1], nodes[4], nodes[3]]]))
    triangles = np.vstack((triangles, [[nodes[3], nodes[4], nodes[6]]]))
    triangles = np.vstack((triangles, [[nodes[3], nodes[6], nodes[7]]]))
    triangles = np.vstack((triangles, [[nodes[7], nodes[6], nodes[5]]]))
    triangles = np.vstack((triangles, [[nodes[7], nodes[5], nodes[0]]]))
    triangles = np.vstack((triangles, [[nodes[5], nodes[6], nodes[4]]]))
    triangles = np.vstack((triangles, [[nodes[5], nodes[4], nodes[2]]]))
    triangles = np.vstack((triangles, [[nodes[3], nodes[7], nodes[0]]]))
    triangles = np.vstack((triangles, [[nodes[3], nodes[0], nodes[1]]]))

    window = MainWindow(root_widget, triangles)

    window.updateTriangles()

    root_widget.mainloop()

main()