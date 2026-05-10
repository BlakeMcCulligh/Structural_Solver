
import scipy.optimize as opt

from OLD_StructuralAnalysis.CrossSectionOptimization.GeneralFunctions import createBounds

def objectiveFunction(X, Constants):
    [solverObject] = Constants
    solverObject.optimizeSolve(X)

    deflection = solverObject.U

    cost = max(deflection)

    return cost

def optimize(solverObject, rangeA, initalGuess):

    Constants = [solverObject]

    bounds = createBounds(rangeA)

    X = initalGuess

    OptimizeResult = opt.minimize(objectiveFunction, X, args=(Constants), method='SLSQP', tol=0.001, bounds=bounds)

    return OptimizeResult.x