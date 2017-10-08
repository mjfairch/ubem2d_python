import unittest
import math
import numpy as np
import numpy.linalg as nla
import ubem2d as ubem

class test_circular_cylinder_flow(unittest.TestCase):
    def test_circular_cylinder_flow(self):
        # Define onset flow and geometry
        uinf = (1,0)
        for n in [4,10,25,50]:
            cyl = ubem.CircularCylinder(n)

            # Solve for stead flow
            soln = ubem.solve_source_body(uinf,cyl)

            # Convenience variables
            sigma = soln.sigma
            cp = soln.cp
            At = soln.At
            An = soln.An
            qn = soln.qn

            # Test shape of solution variables
            self.assertEqual(sigma.shape,(n,))
            self.assertEqual(cp.shape,(n,))
            self.assertEqual(An.shape,(n,n))
            self.assertEqual(At.shape,(n,n))

            # Test Neumann boundary conditions
            self.assertAlmostEqual(nla.norm(qn),0.)

            # Test symmetry and skew symmetry of influence matrices
            self.assertTrue(np.allclose(An.transpose(1,0),An))
            self.assertTrue(np.allclose(At.transpose(1,0),-At))

            # Test diagonal (self-influence) entries of influence matrices.
            # (The exact values of 0.5 and 0.0 were computed by hand)
            for i in range(n):
                self.assertEqual(An[i,i],.5)
                self.assertEqual(At[i,i],0.)

            # Test that source distribution is symmetric and sums to zero
            for i in range(math.floor(n/2)):
                self.assertAlmostEqual(cp[i],cp[-(i+1)])
            self.assertAlmostEqual(np.sum(sigma),0.)

            # Test if pressure agrees with the exact theoretical result
            th = np.arctan2(cyl.ymid, cyl.xmid)
            cp_exact = 1.-4.*np.sin(th)**2
            self.assertAlmostEqual(nla.norm(cp-cp_exact),0.)

            # Test that net forces and moments are zero
            (CD,CL,CM) = ubem.body_cdclcm(uinf, cyl, cp)
            self.assertAlmostEqual(CD,0)
            self.assertAlmostEqual(CL,0)
            self.assertAlmostEqual(CM,0)

if __name__ == '__main__':
    unittest.main()
