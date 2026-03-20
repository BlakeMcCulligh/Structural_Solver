
import scipy.optimize as opt

def createBounds(numCrossSections, rangeA):
    if isinstance(rangeA[0], float):
        b = opt.Bounds([rangeA[0]]*numCrossSections, [rangeA[1]]*numCrossSections)
    else:
        b = opt.Bounds(rangeA[0], rangeA[1])
    return b

def objectiveFunction(X, Constants):
    [solverObject] = Constants
    solverObject.optimizeSolve(X)

    deflection = solverObject.U

    cost = max(deflection)

    return cost

def optimize(solverObject, rangeA, initalGuess):

    Constants = [solverObject]

    bounds = createBounds(len(initalGuess), rangeA)

    X = initalGuess

    OptimizeResult = opt.minimize(objectiveFunction, X, args=(Constants), method='SLSQP', tol=0.001, bounds=bounds)

    return OptimizeResult.x