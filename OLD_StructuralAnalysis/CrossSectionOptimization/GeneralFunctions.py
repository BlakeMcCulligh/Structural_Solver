
import scipy.optimize as opt

def createBounds(rangeVarables):
    return opt.Bounds(rangeVarables[0], rangeVarables[1])