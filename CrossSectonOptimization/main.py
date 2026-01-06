import numpy as np
import math

# cost function - equation 1
def cost(Ai, Lim, row_i, Nd, Nmi):
    """

    :param Ai: the cross-sectional area of members in group i
    :param Lim: the length of member m belonging to group i
    :param row_i: the unit weight of members in group i
    :param Nd: the number of design variables (groups of members with identical cross-sections)
    :param Nmi: the number of members in group i

    :return:
    """

    W = 0
    for i in range(Nd):
        Li = 0
        for m in range(Nmi):
            Li += Lim[i][m]

        W += row_i*Ai*Li

    return W

# displacement and fabricational constraints - equation 2 and 3
def displacementAndFabricationalConstraints(K, N, Ns, L, AiL, AiU, ujk, rjL, rjU, rsL, rsU, usk, Nd, Ai):
    """

    :param K: the loading number
    :param N: the total number of displacement degrees of freedom
    :param Ns: the number of interstory drift constraints
    :param L: the number of loading conditions
    :param AiL:  the lower bound on the cross-sectional area
    :param AiU:  the upper bound on cross-sectional area
    :param ujk: the displacement of the j'th degree of freedom due to loading condition k
    :param rjL: the lower bound on the displacement of the j'th degree of freedom
    :param rjU: the upper bound on the displacement of the j'th degree of freedom
    :param rsL: the lower bounds on the interstory drift
    :param rsU: the upper bounds on the interstory drift
    :param usk: the maximum interstory drift at the s'th floor
    :param Nd: the number of design variables (groups of members with identical cross-sections)
    :param Ai: the cross-sectional area of members in group i
    :return:
    """

    inConstraints = True

    for j in range(N):
        for k in range(L):
            if not (rjL[j] <= ujk[j][k] <= rjU[j]):
                inConstraints = False

    for s in range(Ns):
        for k in range(L):
            if not (rsL[s] <= ujk[s][k] <= rsL[s]):
                inConstraints = False

    for i in range(Nd):
        if not (AiL[i] <= Ai[i] <= AiU[i]):
            inConstraints = False

    return inConstraints

# Stress and buckling constraints based on the AISC ASD specifications - equation 4
def StressBucklingConstraints_AISC_ASD(sigma_mk, Fy, Nm, Fam):
    """

    :param sigma_mk: the stress in member m due to the loading condition k
    :param Fy:  the yield stress of the steel
    :param Nm: the total number of members in the structure
    :param Fam: the allowable axial compressive stress given as a function of the slenderness ratio
    :return:
    """

    stressBucklingConstraints = True
    for m in Nm:
        if -Fam[m] <= sigma_mk[m] <= 0.6*Fy:
            stressBucklingConstraints = False
    return stressBucklingConstraints


# equation 5
def TotalNumberOfMembers(Nd, Nmi):
    """

    :param Nd: the number of design variables (groups of members with identical cross-sections)
    :param Nmi: the number of members in group i
    :return: Nm: the total number of members in the structure
    """
    Nm = 0
    for i in range(Nd):
        Nm += Nmi[i]

    return Nm

# equation 6
def allowableAxialCompressiveStress(Lm, rm, Cc, Fy, E, Nm):
    """

    :param Lm: the length of member m
    :param rm: the radius of Gyration
    :param Cc: the critical slenderness ratio that separates elastic from inelastic buckling
    :param Fy: the yield stress of the steel
    :param E: The elastic modulus of steel
    :param Nm: the total number of members in the structure
    :return: Fam: the allowable axial compressive stress given as a function of the slenderness ratio
    """

    Fam = []
    for m in range(Nm):
        if Lm/rm <= Cc:
            Fam.append(((1-(Lm/rm)**2/(Cc**2))*Fy)/(5/3+(3*(Lm/rm))/(8*Cc)-(Lm/rm)**3/(Cc**3)))
        else:
            Fam.append((12*math.pi*E)/(23*(Lm/rm)**3))

    return Fam

def E14(Amp, sigma_mk, Fam, Fy):
    """

    :param Amp:
    :param sigma_mk: the stress in member m due to the loading condition k
    :param Fam: the allowable axial compressive stress given as a function of the slenderness ratio
    :param Fy: the yield stress of the steel
    :return:
    """

    if abs(sigma_mk) == sigma_mk:
        Amp1 = Amp*(sigma_mk/0.6*Fy)
    else:
        Amp1 = Amp*(sigma_mk/Fam)

    return Amp1

def StressBucklingConstraints_AISC_LRFD(Fcm, sigma_mk, Fy, Nm):
    for m in range(Nm):
        if not (-Fcm[m] <= sigma_mk[m] <= 0.6*Fy):
            return False
    return True


def  criticalCompressiveStrength(lambda_cm, Fy, Lm, rm, E, Nm):
    Fcm = []
    for m in range(Nm):
        if lambda_cm[m] <= 1.5:
            Fcm.append((0.658**(lambda_cm[m]**2))*Fy)
        else:
            Fcm.append((0.877/lambda_cm[m]**2)*Fy)

def lambda_cm(Lm, rm, Fy, E, Nm):

    Lambda_cm = []
    for m in range(Nm):
        Lambda_cm.append(Lm[m]/(rm[m]*math.pi)*(Fy/E)**0.5)


    return Lambda_cm

def ifDiplacementConstraintCritical():

    for i in range(Nd):
        Aip[i]*()

def legrangeMultiplyer(Nac, ujk, W):

    lambda_jk = 0
    for c in range(Nac):
        lambda_jk += W/ujk

    gijk =

    -v_imj










