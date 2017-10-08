import unittest
import math
import numpy as np
import random as rd
from ubem2d.motion.SE2 import SE2 

class test_SE2(unittest.TestCase):
    def random(self):
        return SE2(rd.random(),rd.random(),rd.random())
        
    def test_identity_element(self):
        g = SE2.id()
        self.assertTrue(g.theta == 0. and g.x == 0 and g.y == 0)
    
    def test_multiply(self):
        self.assertTrue(SE2(np.pi/4)*SE2(-np.pi/4) == SE2.id())
        self.assertTrue(SE2(0,1,-2)*SE2(0,-1,2) == SE2.id())
        tr = SE2(0,1,0)
        rot = SE2(np.pi/4,0,0)
        self.assertEqual(tr*rot,SE2(np.pi/4,1,0))
        g = rot*tr
        self.assertEqual(g.theta,np.pi/4)
        self.assertAlmostEqual(g.x,.5*np.sqrt(2))
        self.assertAlmostEqual(g.y,.5*np.sqrt(2))
    
    def test_inversion(self):
        for i in range(10):
            g = self.random()
            ginvg = g.inv()*g
            gginv = g*g.inv()
            self.assertAlmostEqual(ginvg.theta,0.)
            self.assertAlmostEqual(ginvg.x,0.)
            self.assertAlmostEqual(ginvg.y,0.)
            self.assertAlmostEqual(gginv.theta,0.)
            self.assertAlmostEqual(gginv.x,0.)
            self.assertAlmostEqual(gginv.y,0.)
    
    def test_cth_sth(self):
        for i in range(10):
            g = self.random()
            self.assertEqual(math.cos(g.theta),g.cth)
            self.assertEqual(math.sin(g.theta),g.sth)
            new_theta = rd.random()
            g.theta = new_theta
            self.assertEqual(math.cos(g.theta),g.cth)
            self.assertEqual(math.sin(g.theta),g.sth)
    
    def test_mat(self):
        for i in range(10):
            g = self.random()
            A = g.mat()
            self.assertEqual(A.shape,(3,3))
            self.assertEqual(A[0,0],g.cth)
            self.assertEqual(A[0,1],-g.sth)
            self.assertEqual(A[0,2],g.x)
            self.assertEqual(A[1,0],g.sth)
            self.assertEqual(A[1,1],g.cth)
            self.assertEqual(A[1,2],g.y)
            self.assertEqual(A[2,0],0.)
            self.assertEqual(A[2,1],0.)
            self.assertEqual(A[2,2],1.)
    
    def test_r3(self):
        for i in range(10):
            g = self.random()
            (theta,x,y) = g.r3()
            self.assertEqual(theta,g.theta)
            self.assertEqual(x,g.x)
            self.assertEqual(y,g.y)
    
    def test_map_vector(self):
        g = SE2.id()
        x,y = rd.random(), rd.random()
        xx,yy = g.map_vector(x,y)
        self.assertEqual(x,xx)
        self.assertEqual(y,yy)

        g = self.random()
        n = 5 
        x = np.random.rand(n)
        y = np.random.rand(n)
        (xx,yy) = g.map_vector(x,y)
        for i in range(n):
            self.assertEqual(xx[i], g.cth*x[i] - g.sth*y[i])
            self.assertEqual(yy[i], g.sth*x[i] + g.cth*y[i])
    
    def test_map_point(self):
        g = SE2.id()
        x,y = rd.random(), rd.random()
        xx,yy = g.map_vector(x,y)
        self.assertEqual(x,xx)
        self.assertEqual(y,yy)

        g = SE2(np.pi/2,0,0)
        x,y = 1.,0.
        xx,yy = g.map_point(x,y)
        self.assertAlmostEqual(xx,0)
        self.assertAlmostEqual(yy,1)

        xx,yy = g.map_point(x,y,x,y)
        self.assertAlmostEqual(xx,x)
        self.assertAlmostEqual(yy,y)

        g = self.random()
        n = 5 
        x0 = rd.random()
        y0 = rd.random()
        x = np.random.rand(n)
        y = np.random.rand(n)
        (xx,yy) = g.map_point(x,y,x0,y0)
        for i in range(n):
            self.assertEqual(xx[i], g.cth*(x[i]-x0)-g.sth*(y[i]-y0)+g.x+x0)
            self.assertEqual(yy[i], g.sth*(x[i]-x0)+g.cth*(y[i]-y0)+g.y+y0)
    
    def test_equality(self):
        for i in range(10):
            g = self.random()
            self.assertTrue(g == g)
    
    def test_addition_commutative(self):
        for i in range(10):
            g1 = SE2(0,np.random.rand(),np.random.rand())
            g2 = SE2(0,np.random.rand(),np.random.rand())
            h1 = g1*g2
            h2 = g2*g1
            self.assertAlmostEqual(h1.theta, h2.theta)
            self.assertAlmostEqual(h1.x, h2.x)
            self.assertAlmostEqual(h1.y, h2.y)
    
    def test_rotation_commutative(self):
        for i in range(10):
            g1 = SE2(np.random.rand(),0,0)
            g2 = SE2(np.random.rand(),0,0)
            h1 = g1*g2
            h2 = g2*g1
            self.assertAlmostEqual(h1.theta, h2.theta)
            self.assertAlmostEqual(h1.x, h2.x)
            self.assertAlmostEqual(h1.y, h2.y)
    
    def test_associative(self):
        for i in range(10):
            g1 = self.random()
            g2 = self.random()
            g3 = self.random()
            h1 = g1*(g2*g3)
            h2 = (g1*g2)*g3
            self.assertAlmostEqual(h1.theta, h2.theta)
            self.assertAlmostEqual(h1.x, h2.x)
            self.assertAlmostEqual(h1.y, h2.y)
    
if __name__ == '__main__':
    unittest.main()
