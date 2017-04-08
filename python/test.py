import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as Spline


# read data
data = np.genfromtxt('../data/data.csv', delimiter=',')
xData = data[:,0]
yData = data[:,1]
weights = data[:,2]
sl = data[:,3]

# smooth parameter
p = 0.003355272021954153

fit = Spline(xData,yData, s = 10000000, w = weights)
pr = fit(xData)

plt.figure(1)
plt.subplot(211)
plt.plot(xData,yData,'ro')
plt.plot(xData,sl)
plt.plot(xData,pr)
plt.xlabel('x')
plt.ylabel('y')
plt.legend(['data','Matlab','Scipy'])

plt.subplot(212)
plt.plot(xData,sl-pr)
plt.legend(['Matlab - Scipy'])
plt.show()