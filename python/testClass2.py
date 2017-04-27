import matplotlib.pyplot as plt
from SmoothBSpline import SmoothBSpline
import numpy as np

if __name__ == '__main__':

    arr = np.loadtxt('../data/data.csv',delimiter = ',')
    print arr
    x = arr[:,0]
    y = arr[:,1]
    w = arr[:,2]
    ym = arr[:,3]


    sp1 = SmoothBSpline(x, y, w, 0)
    sp1.bspl(x,y,w,0.003355272021954153)
    y1 = sp1.eval(x)


    plt.plot(x,y,'go')
    plt.plot(x,ym,'-r')
    plt.show()




