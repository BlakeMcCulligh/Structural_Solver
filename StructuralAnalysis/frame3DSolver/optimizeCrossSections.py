from StructuralAnalysis.frame3DSolver.__main__ import Frame3D

from StructuralAnalysis.CrossSectionCalculaters import Angle, RectHSS, SquareHSS, TubeHSS


import scipy.optimize as opt

def cost(D, DX, DY, DZ, RX, RY, RZ, weight, reactions, internalForces):
    return max(DZ[0]) + sum(weight) **5

def get_cost(X, constants):
    """
    Gets the cost for the current optimization variables.

    :param X: list. Optimization variables.
    :param constants: list. arguments needed to calculate the cost: [frame, memberGroup, memberGroupType, getWeight, getReactions, getInternalForces, log]
    :return: float. cost
    """

    frame, memberGroup, memberGroupType, getWeight, getReactions, getInternalForces, log = constants

    crossSectionProps = getCrossSectionProps(X, memberGroup, memberGroupType)
    if log: print("Variable cross section properties: ", crossSectionProps)

    j = 0
    for i in range(len(frame.members_CrossSectionProps)):
        if not frame.members[i][3]:
            frame.members_CrossSectionProps[i] = crossSectionProps[j]
            j += 1
    if log: print("Cross Seciton Properites: ", frame.members_CrossSectionProps)

    D, DX, DY, DZ, RX, RY, RZ, weight, reactions, internalForces = frame.analysis_linear(getWeight, getReactions, getInternalForces, log)

    return cost(D, DX, DY, DZ, RX, RY, RZ, weight, reactions, internalForces)

def globalOptimization(frame: Frame3D, memberGroup: list, memberGroupType: list, lowerBound: list or float, upperBound: list or float, getWeight = False, getReactions = False, getInternalForces = False, log=False):
    """
    Finds the optimum cross-sections for members with their cross-sections not set.

    :param frame: Frame3D. 3D Frame to be optimized.
    :param memberGroup: list. list of indices of member groups for non set members to be assigned to. Must be length of non set members.
    :param memberGroupType: list. List of cross-section types for each member group to be assigned. Must be length of number of member gorups.
    :param lowerBound: list or float. Lower bound on the optimization variables. if list must be length of number of variables.
    :param upperBound: list or float. Upper bound on the optimization variables. if list must be length of number of variables.
    :param getWeight: bool. Weather the weighht of all the members is needed for the cost function.
    :param getReactions: bool. Weather the reactions are needed for the cost function.
    :param getInternalForces: bool. Weather the internal forces are needed for the cost function.
    :return: scipi optimization_results class: results of the optimization.
    """

    if log:
        print("--------------------------------------------------------------")
        print("-----------------  Global Optimization  ----------------------")
        print("--------------------------------------------------------------")

    chackInputs(frame, memberGroup, memberGroupType)

    numVarables = getNumVarables(memberGroupType)
    if log: print("numVarables: ", numVarables)

    frame.preAnalysis_linear(log=log)

    constants = [frame, memberGroup, memberGroupType, getWeight, getReactions, getInternalForces, log]

    bounds = getBounds(lowerBound, upperBound, numVarables)

    optimization_results = opt.shgo(get_cost, bounds, args=[constants])

    if log: print("optimization results: ", optimization_results)

    return optimization_results


def chackInputs(frame: Frame3D, memberGroup: list, memberGroupType: list):
    """
    Checks if the input lists for the optimization are valid.

    :param frame: Frame3D. 3D Frame to be optimized.
    :param memberGroup: list. list of indices of member groups for non set members to be assigned to.
    :param memberGroupType: list. List of cross-section types for each member group to be assigned.
    """

    numberNotSetMembers = 0
    for member in frame.members:
        if not member[3]:
            numberNotSetMembers += 1
    if len(memberGroup) != numberNotSetMembers:
        raise Exception('Number of not set members is not equal to number of members assigned to a group.')
    if max(memberGroup) + 1 != len(memberGroupType):
        raise Exception('Number of member groups does not equal number of groups members are assigned to')

def getNumVarables(memberGroupType: list):
    """
    Gets the number of variables to be in the optimization problem.

    :param memberGroupType: list. List of cross-section types for each member group to be assigned.
    :return:
        numVariables: int, number of variables in the optimization.
    """

    numVarables = 0
    for t in memberGroupType:
        if t == "Angle":
            numVarables += 3
        elif t == "RectHSS":
            numVarables += 3
        elif t == "SquareHSS":
            numVarables += 2
        elif t == "TubeHSS":
            numVarables += 2
        else:
            raise Exception('member groupe type is not a valid cross section type. Chose one of the following: Angle, RectHSS, SquareHSS, TubeHSS')
    return numVarables

def getBounds(lowerBound: list or float, upperBound: list or float, numVarables: int):
    """
    Sets up the scipy optimization bounds.

    :param lowerBound: list or float. Lower bound on the optimization variables. if list must be length of number of variables.
    :param upperBound: list or float. Upper bound on the optimization variables. if list must be length of number of variables.
    :param numVarables: int, Number of variables to be in the optimization problem.
    :return: Scipy optimization bounds.
    """

    if isinstance(lowerBound, float):
        lowerBound = [lowerBound] * numVarables
    if isinstance (upperBound, float):
        upperBound = [upperBound] * numVarables
    return opt.Bounds(lowerBound, upperBound)

def getCrossSectionProps(X, memberGroups, memberGroupType):
    """
    Gets the cross-section properties of the cross-section being optimized from the optimization variables.

    :param X: list: Optimization variables.
    :param memberGroups: list. list of indices of member groups for non set members to be assigned to.
    :param memberGroupType: list. List of cross-section types for each member group to be assigned.
    :return: list: Cross-section properties for the optimization cross-sections. shape: (# cross-sections being optimized, 4
    """

    groupProperties = []
    j = 0
    for t in memberGroupType:
        if t == "Angle":
            A = Angle.getA(X[j], X[j+1], X[j+2])
            Iy = Angle.getIy(X[j], X[j+1], X[j+2])
            Iz = Angle.getIx(X[j], X[j+1], X[j+2])
            J = Angle.getJ(X[j], X[j+1], X[j+2])
            groupProperties.append([A, Iy, Iz, J])
            j += 3
        elif t == "RectHSS":
            A = RectHSS.getA(X[j], X[j+1], X[j+2])
            Iy = RectHSS.getIy(X[j], X[j+1], X[j+2])
            Iz = RectHSS.getIx(X[j], X[j+1], X[j+2])
            J = RectHSS.getJ(X[j], X[j+1], X[j+2])
            groupProperties.append([A, Iy, Iz, J])
            j += 3
        elif t == "SquareHSS":
            A = SquareHSS.getA(X[j], X[j+1])
            I = SquareHSS.getI(X[j], X[j+1])
            J = SquareHSS.getJ(X[j], X[j+1])
            groupProperties.append([A, I, I, J])
            j += 2
        elif t == "TubeHSS":
            A = TubeHSS.getA(X[j], X[j+1])
            I = TubeHSS.getI(X[j], X[j+1])
            J = TubeHSS.getJ(X[j], X[j+1])
            groupProperties.append([A, I, I, J])
            j += 2

    memberProperties = []
    for memberG in memberGroups:
        memberProperties.append(groupProperties[memberG])

    return memberProperties