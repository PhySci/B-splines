import numpy as np
import matplotlib.pyplot as plt

def sortAvg(x,y,w):
    """
    Sorting and average of degenerated knots
    :param x: abscissa array 
    :param y: ordinata array
    :param w: weight array
    :return: list of new x,y,w arrays
    """

    z = np.diff(x)
    z[z != 0] = 1
    z = np.cumsum(np.hstack([0, z]))

    maxInd = int(z[-1])+1
    xNew = np.zeros(maxInd)
    yNew = np.zeros(maxInd)
    wNew = np.zeros(maxInd)

    for ind in range(maxInd):
        degInd = np.argwhere(z == ind)
        if (degInd.shape[0] == 1):
            xNew[ind] = x[degInd]
            wNew[ind] = w[degInd]
            yNew[ind] = y[degInd]
        else:
            xNew[ind] = x[degInd[0]]
            wNew[ind] = np.mean(w[degInd])
            yNew[ind] = np.sum(y[degInd]*w[degInd])/np.sum(w[degInd])

    return xNew, yNew, wNew


if __name__=='__main__':
    nPoints = 70
    x = np.linspace(0,10,nPoints)

    y = np.cos(x)
    # introduce multiplicity of knots
    x[5] = x[6]
    x[30] = x[31]
    x[60] = x[59]

    w = np.ones(nPoints)
    w[5] = 0.4
    w[6] = 0.6
    w[30] = 0.8
    w[31] = 1.2
    w[60] = 1.5
    w[59] = 0.8


    [xNew,yNew,wNew] = sortAvg(x,y,w)

    plt.figure(1)
    plt.subplot(211)
    plt.plot(x,y,'ro')
    plt.plot(xNew,yNew,'gx')

    plt.subplot(212)
    plt.plot(x,w,'ro')
    plt.plot(xNew,wNew,'gx')
    plt.show()
