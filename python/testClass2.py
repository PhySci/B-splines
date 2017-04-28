import matplotlib.pyplot as plt
from SmoothBSpline import SmoothBSpline
import numpy as np

if __name__ == '__main__':

    arr = np.loadtxt('../data/data.csv',delimiter = ',')


    x = arr[:,0]
    y = arr[:,1]
    w = arr[:,2]
    ym = arr[:,3]



    sp1 = SmoothBSpline()
    p = 0.003355272021954153
    sp1.bspl(x,y,w,p)
    y1 = sp1.eval(x)

    sp2 = SmoothBSpline()
    p = 0.003355272021954153
    sp2.bspl(x,y,None,p)
    y2 = sp2.eval(x)

    plt.figure(1)
    plt.subplot(211)
    plt.plot(x,y,'go')
    plt.plot(x,ym,'-r')
    plt.plot(x,y1,'--k')
    plt.legend(['experimental points','MatLab','Python'])
    plt.title(['p =',p])

    plt.subplot(212)
    plt.semilogy(x,np.abs(ym-y1))
    plt.legend(['MatLab - Python'])

    plt.figure(2)
    plt.subplot(211)
    plt.plot(x, y, 'go')
    plt.plot(x, ym, '-r')
    plt.plot(x, y2, '--k')
    plt.legend(['experimental points', 'MatLab', 'Python'])
    plt.title(['p =', p])

    plt.subplot(212)
    plt.semilogy(x, np.abs(ym - y1))
    plt.legend(['MatLab - Python'])

    plt.figure(3)
    plt.subplot(211)
    plt.plot(x, y1, '-r')
    plt.plot(x, y2, '-g')
    plt.legend(['With weights', 'Without weights'])
    plt.title(['p =', p])

    plt.subplot(212)
    plt.plot(x, y1-y2)


    plt.show()




