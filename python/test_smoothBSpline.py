import unittest
from SmoothBSpline import SmoothBSpline
import numpy as np
import numpy.testing as npt


class TestSmoothBSpline(unittest.TestCase):

    def test_CheckDataEqual(self):
        spl = SmoothBSpline();

        n = 100;
        x = np.linspace(0,n)
        y = np.cos(x)
        w = np.ones(n)
        [a,b,c] = spl.checkData(x, y, w, 0.1)
        self.assertEqual(a, x, 'Not equal')
        self.assertEqual(b, y, 'Not equal')
        self.assertEqual(c, w, 'Not equal')

    def test_checkData(self):

        spl = SmoothBSpline()
        self.assertRaises(Exception,spl.checkData,1,2,3)



    #def test_demo(self):

    #    spl = SmoothBSpline()
    #    self.fail()


if __name__ == '__main__':
    unittest.main()
