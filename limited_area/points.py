from __future__ import absolute_import, division, print_function
import os
import sys

""" Parse the points file, and output a representation that can be passed back
that can be used to.

File Examples:

# Circle
```
Name: Colorado
Type: Circle
center: 38.9386, -105.6875
radius: 385.5 km
```

# Square
```
Name: CONUS
Type: Square
left_lon: -126.0
right_lon: -64.0
up_lat: 50.0
lower_lat: 25.0
```

# Custom
```
Name: pnw
Type: Custom
Point: 45.0, -118.0
50.000, -129.000
50.000, -115.000
41.500, -115.0000
41.500, -129.000
```

# Ellipse
```
Name: Russia
Type: Ellipse
```

# Example
```
Name: 
Type: 
Point: Lat, Lon
Lat1, Lon1
Lat2, Lon2
Lat3, Lon3
Lat4, Lon4
...
LatN, LonN
"""

def PointsParser(self, file, *args, **kwargs):
    """ Parse file for our points syntax and set variables 
    accordingly (self.variable ...)"""

    self.points = []

    # Kwargs
    self._DEBUG_ = kwargs.get('DEBUG', 0)

    if not os.path.isfile(file):
        print("ERROR: This points file could not be found")
        sys.exit(-1)
    else:
        self.points_file = open(file, 'r')

    line_number = 0
    for line in self.points_file:
        line_number += 1

        if '#' in line:
            line = line.split('#')[0]
            line = line.strip()

        if ':' in line:
            lhs = line.split(':')[0]
            rhs = line.split(':')[1]
            lhs = lhs.strip()
            rhs = rhs.strip()

            if lhs == 'Name' or lhs == 'name':
                # TODO: Do some error checking of the name
                self.name = rhs
            # Check to see if the type matches our avaiable types
            elif lhs == 'Type' or lhs == 'type':
                if rhs == 'Custom' or rhs == 'custom':
                    self.type = 'custom'
                elif rhs == 'Square' or rhs == 'square':
                    self.type = 'square'
                elif rhs == 'Circle' or rhs == 'circle':
                    self.type = 'circle'
                elif rhs == 'Ellipse' or rhs == 'ellipse':
                    self.type = 'ellipse'
                else:
                    print("ERROR: This is not a valid points type: ", rhs)
                    sys.exit(-1)
            elif lhs == 'Point' or lhs == 'point':
                if ',' in rhs:

                    self.in_point = [float(rhs.split(',')[0]),
                                     float(rhs.split(',')[1])]
            elif lhs == 'Radius' or lhs == 'radius':
                self.radius = float(rhs)

        elif line == 'keyword': # Then we have a keyword option
            pass
        elif ',' in line: # Then we have a coordinate point
            lat = float(line.split(',')[0])
            lon = float(line.split(',')[1])

            self.points.append(lat)
            self.points.append(lon)

    self.points_file.close()
    del(self.points_file)
       

if __name__ == "__main__":
    print(sys.argv[1])
    parse_points(sys.argv[1])

