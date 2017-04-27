import matplotlib.pyplot as plt
from SmoothBSpline import SmoothBSpline
import numpy as np

if __name__ == '__main__':

    nPoints = 100
    x = - np.linspace(0,10,nPoints)
    y = np.cos(x)
    w = np.ones(nPoints)

    sp1 = SmoothBSpline(x, y, w, 0)
    [x,y,w] = sp1.checkData(x,y,w)

    sp1.bspl(x,y,w,0.1)
    y1 = sp1.eval(x)


    plt.figure(1)
    plt.plot(x,y,'go')
    plt.plot(x,y1,'-r')
    plt.show()