from __future__ import absolute_import, division, print_function

import os
import sys
import types
import numpy as np

from limited_area.shape_reader import ShapeReader
from limited_area.points import PointsParser
import numpy as np

NOT_IMPLEMENTED_ERROR = "IS NOT IMPLEMENTED - YOU SHOULD IMPLENTED IT!"

""" region_spec.py - Provide a number of operations for defining a 
region. """

if sys.version_info[0] > 2:
    create_bound_method = types.MethodType
else:
    def create_bound_method(func, obj):
        return types.MethodType(func, obj, obj.__class__)



def normalize_cords(lat, lon):
    """ Returned lat, and lon to be in radians and the same 
    range as MPAS - Lat: -pi/2 to pi/2 - LonL 0 to 2*pi
       
    Lat - Latitude in degrees
    Lon - Lontitude in degrees

    """
    lat *= np.pi / 180.0
    lon *= np.pi / 180.0

    if lon < 0:
        lon += 2 * np.pi

    return lat, lon


class RegionSpec:
    """ RegionSpec - Method for generating regional specifications
    Region spec works upon a contract. It will need the following 
    information when it is called:

    filename  - Path to the file that is to be read 
    file type - the type of the file that it is to be read.
                Currently aviable options are:
                    * points     - Implemented
                    * Shape File - No Implemented

    And will then return, the following:

    filename         - Output filename (if desired and specified in the specification file)
    points array     - A 1-dimensional list of lat lon cords specificying boundry 
                       in counter clockwise order
    in-point         - pair of points that inside the boundary
    algorithm choice - The desired algorithm for choosing bondary points and 
                       relaxation layers (if specified within the specification file)
    """
    # TODO: Update __init__ with fileName
    def __init__(self, method='points', *args, **kwargs):
        """ init for region Spec

        From the choosen method. Bind that method to this object,
        that way, it will act as our own method.

        method  -- Method used for specifying region. Options are:
                    * points or point
                    * shapeFile or shape

        Keyword Arguments: 
            DEBUG - Debug value for verbose output - Default 0
        """
        # Keyword Args
        self._DEBUG_ = kwargs.get('DEBUG', 0)

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

    def gen_spec(self, fileName, *args, **kwargs):
        """ Generate the specifications and return, name, in point and a list of points.

        Call the method we bound above, and then do any proccessing here
        to do things like convert cordinates to radians, or anything else
        we need to get return contract variables.
        
        fileName - The file that specifies the region

        Return values:
            name     - The name of the region
            in_point - A point that is within the region
            points   - A 1-Dimensional list of latitude and longitude points
                       (in degrees) that list the boundary points of the region
                       in counter-clockwise. ie: [lat1, lon1, lat2, lon2, ... , latN, lonN]
        """

        self._gen_spec(fileName, *args, **kwargs)

        if self.method == 'SHAPE':
            # Do shape stuff here and return the following
            return self.name, self.in_point, self.points

        if self.method == 'POINTS':
            if self.type == 'custom':
                self.points = np.array(self.points)

                # Convert the points to radians and set them to be between
                # Lon: 0 to 2*Pi and Lat: -pi to +pi
                for cord in range(0, len(self.points), 2):
                     self.points[cord], self.points[cord+1] = normalize_cords(
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
        """ """
        raise NotImplementedError("CIRCLE FUCTION "+NOT_IMPLEMENTED_ERROR)

    def square(self):
        """ """
        raise NotImplementedError("SQAURE FUNCITON "+NOT_IMPLEMENTED_ERROR)

    def ellipse(self):
        """ """
        raise NotImplementedError("ELLIPSE FUNCITON "+NOT_IMPLEMENTED_ERROR)


