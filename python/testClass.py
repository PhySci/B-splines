import matplotlib.pyplot as plt
from SmoothBSpline import SmoothBSpline
import numpy as np

if __name__ == '__main__':


    x = np.linspace(0,5,10)
    y = x**3+1
    w = np.ones(10)
    #print 'w=',w
    #print 'y=',y

    spl = SmoothBSpline(x, y, w, 0)
    spl.bspl(x,y,w,0.1)
    #x2 = np.linspace(1,4.99,100)
    #k = spl.calcCurve(x2)


    #plt.plot(x, y,'rx')
    #plt.plot(x2, k, '-g')

    #plt.show()

