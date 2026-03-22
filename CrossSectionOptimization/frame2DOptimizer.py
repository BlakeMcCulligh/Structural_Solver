
import scipy.optimize as opt

from CrossSectionOptimization.GeneralFunctions import createBounds

def objectiveFunction(X, Constants):
    [solverObject] = Constants
    solverObject.optimizeSolve(X)

    deflection = solverObject.U

    cost = max(deflection)

    return cost

def optimize(solverObject, rangeVarables, initalGuess):

    Constants = [solverObject]

    bounds = createBounds(rangeVarables)

    X = initalGuess

    OptimizeResult = opt.minimize(objectiveFunction, X, args=(Constants), method='SLSQP', tol=0.001, bounds=bounds)

    return OptimizeResult.x