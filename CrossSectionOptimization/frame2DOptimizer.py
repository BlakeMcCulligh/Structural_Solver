
import scipy.optimize as opt

def createBounds(rangeVarables):
    return opt.Bounds(rangeVarables[0], rangeVarables[1])

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