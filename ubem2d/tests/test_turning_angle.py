import unittest
import numpy as np
from ubem2d.math.turning_angle import turning_angle

class test_turning_angle(unittest.TestCase):
    def test_ccw_square_open(self):
        x = np.array([0.,1,1,0,0])
        y = np.array([0.,0,1,1,0])
        self.assertEqual(turning_angle(x,y),3*np.pi/2)

    def test_ccw_square_closed(self):
        x = np.array([0.,1,1,0,0])
        y = np.array([0.,0,1,1,0])
        self.assertEqual(turning_angle(x,y,True),2*np.pi)

    def test_cw_square_open(self):
        x = np.array([0.,0,1,1,0])
        y = np.array([0.,1,1,0,0])
        self.assertEqual(turning_angle(x,y),-3*np.pi/2)

    def test_cw_square_closed(self):
        x = np.array([0.,0,1,1,0])
        y = np.array([0.,1,1,0,0])
        self.assertEqual(turning_angle(x,y,True),-2*np.pi)

    def test_ccw_circle_closed(self):
        t = np.linspace(0,2*np.pi,100)
        # Starting at (1,0)
        x = np.cos(t)
        y = np.sin(t)
        self.assertAlmostEqual(turning_angle(x,y,True), 2*np.pi)
        # Starting at (0,1)
        x = -np.sin(t)
        y = np.cos(t)
        self.assertAlmostEqual(turning_angle(x,y,True), 2*np.pi)
        # Starting along branch cut at (-1,0)
        x = -np.cos(t)
        y = -np.sin(t)
        self.assertAlmostEqual(turning_angle(x,y,True), 2*np.pi)
        # Starting at (0,-1)
        x = np.sin(t)
        y = -np.cos(t)
        self.assertAlmostEqual(turning_angle(x,y,True), 2*np.pi)

    def test_cw_circle_closed(self):
        t = np.linspace(0,-2*np.pi,100)
        # Starting at (1,0)
        x = np.cos(t)
        y = np.sin(t)
        self.assertAlmostEqual(turning_angle(x,y,True), -2*np.pi)
        # Starting at (0,1)
        x = -np.sin(t)
        y = np.cos(t)
        self.assertAlmostEqual(turning_angle(x,y,True), -2*np.pi)
        # Starting along branch cut at (-1,0)
        x = -np.cos(t)
        y = -np.sin(t)
        self.assertAlmostEqual(turning_angle(x,y,True), -2*np.pi)
        # Starting at (0,-1)
        x = np.sin(t)
        y = -np.cos(t)
        self.assertAlmostEqual(turning_angle(x,y,True), -2*np.pi)
    
    def test_L_ccw(self):
        # Test the following L-shaped polygon:
        #   ___ 
        #   | |
        # __| |
        # |   |
        # -----
        x = np.array([0.,2,2,1,1,0,0])
        y = np.array([0.,0,2,2,1,1,0])
        self.assertAlmostEqual(turning_angle(x,y,True), 2*np.pi)
    
    def test_random_closed_is_multiple_2pi(self):
        # Take a random collection of points, form the broken line
        # joining them, and then join the last point to the first point, thus
        # forming a closed broken line.  The closed turning angle of the 
        # resulting closed broken line should be an integer multiple of 2*pi.
        n = 1000
        x = np.random.rand(n)
        y = np.random.rand(n)
        x = np.concatenate([x,[x[0]]])
        y = np.concatenate([y,[y[0]]])
        a1 = turning_angle(x,y,True)/(2*np.pi)
        self.assertAlmostEqual(a1-round(a1),0.)
    
    def test_reverse_minus(self):
        # If points are traversed in the opposite direction, the turning
        # angle ought to flip sign.
        x = np.random.rand(25)
        y = np.random.rand(25)
        a1 = turning_angle(x,y)
        a2 = turning_angle(np.flipud(x),np.flipud(y))
        self.assertAlmostEqual(a1,-a2)
