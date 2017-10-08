from enum import Enum
import numpy as np

__all__ = ['Orientation']

class Orientation(Enum):
    CCW = True  # counterclockwise
    CW = False  # clockwise

    def reverse(self):
        if (self == Orientation.CCW):
            return Orientation.CW
        elif (self == Orientation.CW):
            return Orientation.CCW
        else:
            raise TypeError('Invalid orientation')
    
    def __repr__(self):
        if (self == Orientation.CCW):
            return 'Counterclockwise'
        elif (self == Orientation.CW):
            return 'Clockwise'
        else:
            raise TypeError('Invalid orientation')
