import numpy as np

# 2 Problem Formulation

# rows: cross-section
# columbs: member
X = np.array([[]])

# x i, j = 1 if candidate i is chosen in member j, = 0 Otherwise


# minimise f(m(X), C(x))

#subject to
# (1) [K]D = R
# (2) sigma VM <= sigma max
# (3) np.sum(np_array, axis=0) == 1
# (4) x i,j == 1 or 0



# 3 Methodology

# 3.1 Relaxation and interpolation

# change constraint (4) to x i,j == [0,1]

# interpolation
def SIMP(x, p):
    """

    :param x: design variable: amout a member is chosen as a particular cross-section
    :param p: Penalty parameter
                if p >  1: penalization is below-linear
                if 0 < p < 1: penalization is above-linear. NOT USED due to gradian apoching inf as p aproches 0
    :return:
    """

    return x**p

def RAMP(x, r):
    """

    :param x: design variable: amout a member is chosen as a particular cross-section
    :param r: Penalty parameter
                if r > 0: penalization is below-linear
                if -1 < r < 0: penalization is above-linear
    :return:
    """

    #TODO figure out seq

    return x/(1+r*(1-x))


# The values of the penalty parameters are controlled through a continuation method
#
# The penalty is chosen such that the interpolation is
# initially linear or slightly penalized. Through iterations the
# penalty is increased such that the final design hopefully is
# an integer solution. In each iteration the optimized design
# from the previous iteration is used as a starting point. If the
# penalty starts too high, the risk of ending in a local minimum
# is increased.

def DMO_Element_Stiffness_matix(K, K0, ):
    a=1