import unittest
import math
import numpy as np
import numpy.linalg as nla
import ubem2d as ubem

class test_naca_flow(unittest.TestCase):
    def test_symmetric_naca(self):
        # Define onset flow and geometry
        Uinf = (1,0)
        for n in [25,50]:
            codes = ['0001','0005','0010','0015','0020']
            for code in codes:
                # Create airfoil
                foil = ubem.naca4(code,n)

                # Solve steady-flow problem
                soln = ubem.solve_hess_smith_body(Uinf,foil)
                # Convenience variables
                sigma = soln.sigma
                gamma = soln.gamma
                cp = soln.cp
                qn = soln.qn
                At = soln.At
                An = soln.An
                Bt = soln.Bt
                Bn = soln.Bn

                # Test shape of solution variables
                self.assertEqual(sigma.shape,(n,))
                self.assertTrue(np.isscalar(gamma))
                self.assertEqual(cp.shape,(n,))
                self.assertEqual(An.shape,(n,n))
                self.assertEqual(At.shape,(n,n))
                self.assertEqual(Bn.shape,(n,n))
                self.assertEqual(Bt.shape,(n,n))

                # Test diagonal (self-influence) entries of influence matrices
                # (The exact values were computed by hand)
                for i in range(n):
                    self.assertEqual(An[i,i],.5)
                    self.assertEqual(At[i,i],0.)
                    self.assertEqual(Bn[i,i],0)
                    self.assertEqual(Bt[i,i],.5)

if __name__ == '__main__':
    unittest.main()
