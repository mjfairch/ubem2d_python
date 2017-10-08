import unittest
import math
import ubem2d as ubem

class test_steady_lift_model(unittest.TestCase):
    def test_steady_lift_model_symmetric(self):
        codes = ['0001','0005','0010','0015']
        npan = 50
        for code in codes:
            foil = ubem.naca4(code, npan)
            (CL0,m,aoa,CL) = ubem.steady_lift_model(foil,-10,10,11)
            m_theoretical = math.pi**2/90
            # Test if relative error between computed slope and theoretical
            # slope is small
            self.assertTrue(math.fabs(m-m_theoretical)/m < .2)

    def test_steady_lift_model_cambered(self):
        codes = ['1105','2210','4415','6520']
        npan = 50
        for code in codes:
            foil = ubem.naca4(code,npan)
            (CL0,m,aoa,CL) = ubem.steady_lift_model(foil,-10,10,11)
            m_theoretical = math.pi**2/90
            # Test if relative error between computed slope and theoretical
            # slope is small
            self.assertTrue(math.fabs(m-m_theoretical)/m < .2)
