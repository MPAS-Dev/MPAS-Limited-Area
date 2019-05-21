import os
import sys
import types
import numpy as np

'''
# class region_spec 

```
To Region Spec
==============
filename - Path to the file that is to be read 
file type - the type of the file that it is to be read - How will this be
            specirfied?
```


```
Back To Limited Area 
==================
filename - Output filename (if desired)
points array - 2d list of lat lon cords specificying boundry
in-point - pair of points that inside the boundary
algorithm choice - if any
```


'''

from LimitedArea.shapeReader import ShapeReader
from LimitedArea.points import PointsParser
import numpy as np


if sys.version_info[0] > 2:
    create_bound_method = types.MethodType
else:
    def create_bound_method(func, obj):
        return types.MethodType(func, obj, obj.__class__)


NOT_IMPLEMENTED_ERROR = "IS NOT IMPLEMENTED - YOU SHOULD IMPLENTED IT!"

# Normalize cords to be:
# 1. In radians
# 2. Lat = -pi - pi
# 3. Lon = 0 - 2*pi
def normalize_cords(lat, lon):
    lat *= np.pi / 180.0
    lon *= np.pi / 180.0

    if lon < 0:
        lon += 2 * np.pi

    return lat, lon


class RegionSpec:
    def __init__(self, method='points', *args, **kwargs):

        if method == 'Points' or method == 'points':
            self._gen_spec = create_bound_method(PointsParser, self)
            self.method = 'POINTS'
        elif (method == 'shape' 
          or method == 'shapefile' 
          or method == 'shapeFile'):
            self._gen_spec = create_bound_method(ShapeReader, self)
            self.method = 'SHAPE'
        else:
            raise NotImplementedError("SHAPE READER "+NOT_IMPLEMENTED_ERROR)

    def gen_spec(self, file, *args, **kwargs):
        self._gen_spec(file, *args, **kwargs)

        if self.method == 'SHAPE':
            # Do shape stuff here and return the following
            return self.name, self.in_point, self.points

        if self.method == 'POINTS':
            # Do Points stuff here and return
            if self.type == 'custom':
                self.points = np.array(self.points)

                for cord in range(0, len(self.points), 2):
                     (self.points[cord], self.points[cord+1]) = normalize_cords(
                                                                self.points[cord], 
                                                                self.points[cord+1])

                self.in_point[0], self.in_point[1] = normalize_cords(
                                                        self.in_point[0],
                                                        self.in_point[1])

                return self.name, self.in_point, self.points
            elif self.type == 'square':
                return self.square()
            elif self.type == 'circle':
                return self.circle()
            elif self.type == 'ellipse':
                return self.ellipse()


    def circle(self):
        raise NotImplementedError("CIRCLE FUCTION "+NOT_IMPLEMENTED_ERROR)

    def square(self):
        raise NotImplementedError("SQAURE FUNCITON "+NOT_IMPLEMENTED_ERROR)

    def ellipse(self):
        raise NotImplementedError("ELLIPSE FUNCITON "+NOT_IMPLEMENTED_ERROR)


