import matplotlib.pyplot as plt
from python.SmoothBSpline import SmoothBSpline
import numpy as np

if __name__ == '__main__':


    n = 100
    x = np.linspace(0,10,n)
    y = np.cos(x)
    w = np.ones(n)


    sp1 = SmoothBSpline(x, y, w, 0)
    sp1.bspl(x,y,w,1)
    y1 = sp1.eval(x)

    sp2 = SmoothBSpline(x, y, w, 0)
    sp2.bspl(x,y,w,0.5)
    y2 = sp2.eval(x)

    sp3 = SmoothBSpline(x, y, w, 0)
    sp3.bspl(x, y, w, 0.1)
    y3 = sp3.eval(x)

    sp4 = SmoothBSpline(x, y, w, 0)
    sp4.bspl(x, y, w, 0)
    y4 = sp4.eval(x)

    plt.subplot(311)
    plt.plot(x, y,'rx')
    plt.plot(x, y1, '-b')
    plt.legend(['1'])

    plt.subplot(312)
    plt.plot(x, y,'rx')
    plt.plot(x, y2, '-b')
    plt.legend(['0.5'])

    plt.subplot(313)
    plt.plot(x, y,'rx')
    plt.plot(x, y3, '-b')
    plt.legend(['0.1'])

    plt.show()