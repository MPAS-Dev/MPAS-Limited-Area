from __future__ import absolute_import, division, print_function

import os
import sys
import types
import numpy as np

from limited_area.shape_reader import ShapeReader
from limited_area.points import PointsParser
from limited_area.mesh import latlon_to_xyz, xyz_to_latlon
from limited_area.mesh import rotate_about_vector
import numpy as np

NOT_IMPLEMENTED_ERROR = "IS NOT IMPLEMENTED - YOU SHOULD IMPLENTED IT!"

EARTH_RADIUS = 6371229 # Meters

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
                if self._DEBUG_ > 0:
                    print("DEBUG: Using a custom poloygon for generating a region")

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
                if self._DEBUG_ > 0:
                    print("DEBUG: Using the circle method for region generation")

                self.in_point[0], self.in_point[1] = normalize_cords(
                                                        self.in_point[0],
                                                        self.in_point[1])

                # Convert to meters, then divide by radius to get radius upon sphere w/ r = 1
                self.radius = (self.radius * 1000) / EARTH_RADIUS
                self.points = self.circle(self.in_point[0], self.in_point[1], self.radius)

                return self.name, self.in_point, self.points.flatten()

            elif self.type == 'ellipse':
                return self.ellipse()



    def circle(self, center_lat, center_lon, radius):
        """ Return a list of latitude and longitude points in degrees that
        area radius away from (center_lat, center_lon)

        center_lat - Circle center latitude in radians
        center_lon - Circle center longitude in radians
        radius     - Radius of desire circle in radians upon the unit shpere

        """

        P = []

        if self._DEBUG_ > 1:
            print("DEBUG: center_lat: ", center_lat,
                        " center_lon: ", center_lon,
                        " radius: ", radius)


        C = latlon_to_xyz(center_lat, center_lon, 1.0)

        # Find a point not equal to C or -C
        K = np.zeros(3)
        eps = 0.1

        if ((C[0] >= eps or C[0] <= eps) and (C[1] >= eps or C[0] <= eps)):
            K[0] = 0.0
            K[1] = 1.0
            K[2] = 0.0
        else:
            K[0] = 1.0
            K[1] = 0.0
            K[2] = 0.0

        # S is then a vector orthogonal to C
        S = np.cross(C, K)

        P0 = rotate_about_vector(C, S, radius)

        for r in np.linspace(0.0, 2.0*np.pi, 100):
            P.append(rotate_about_vector(P0, C, r))

        ll = []
        # TODO: The efficency here can be improved for memory
        # and probably comp time
        for i in range(len(P)):
            ll.append(xyz_to_latlon(P[i])) # Convert back to latlon

        return np.array(ll)


    def square(self):
        """ """
        raise NotImplementedError("SQAURE FUNCITON "+NOT_IMPLEMENTED_ERROR)


    def ellipse(self):
        """ """
        raise NotImplementedError("ELLIPSE FUNCITON "+NOT_IMPLEMENTED_ERROR)
