import unittest
from python.SmoothBSpline import SmoothBSpline
import numpy as np
import numpy.testing as npt
import os.path
import matplotlib.pyplot as plt


class TestSmoothBSpline(unittest.TestCase):

    @classmethod
    def test_CheckDataEqual(self):
        spl = SmoothBSpline()
        n = 100
        x = np.linspace(0,10,n)
        y = np.cos(x)
        w = np.ones(n)
        [a,b,c] = spl.checkData(x, y, w, 0.1)
        npt.assert_array_equal(a, x, 'Not equal')
        npt.assert_array_equal(b, y, 'Not equal')
        npt.assert_array_equal(c, w, 'Not equal')

    def test_checkData(self):

        spl = SmoothBSpline()
        #self.assertRaises(Exception,spl.checkData,1,2,3)

    def test_simpleSet(self):
        x = np.linspace(0, 5, 10)
        y = x ** 3
        sp1 = SmoothBSpline()
        sp1.bspl(x, y, None, 1)
        y1 = sp1.eval(x)
        npt.assert_almost_equal(y1, y, -2)


    def test_realDataSet(self):
        path = os.path.join(os.path.dirname(__file__),'data','data.csv')
        arr = np.loadtxt(path, delimiter=',')

        x = arr[:, 0]
        y = arr[:, 1]
        w = arr[:, 2] # weight
        ym = arr[:, 3] # matlab prediction

        sp1 = SmoothBSpline()
        p = 0.003355272021954153
        sp1.bspl(x, y, w, p)
        y1 = sp1.eval(x)
        npt.assert_almost_equal(y1,ym,0)

if __name__ == '__main__':
    unittest.main()
