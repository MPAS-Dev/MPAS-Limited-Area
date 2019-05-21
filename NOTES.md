
# MGDUDA Notes & Algorithms

**High Level**
1. Mark Cells in Region - See below
2. Build relaxtion layers - Easy similar to 3
3. Build specified layers - Easy similar to 2
4. Mark all edges and verties in region - also easy
5. Reindexing - Renumber cells, edges, verties to a continuous range
    * Perhaps a little harder
    * Build a mapping from global ids to region ids
    * Then remape cells on cells and edges on edge
    * Reindexing
    * Perhaps generic for fields
6. Loop over input fields, writing subsets to output file

**Marking Cells in Region**
1.1. Get list of points approximating the boundary
1.2. Mark paths connecting these
1.3. Floodfill - Need point inside 

* Read in sphere units from mesh

# Classes

1. Limited Area
2. Region Spec
3. Points Parser - Maybe a subclass of region spec?
4. Shape Reader - Maybe a subclass of region spec?

interface between limited area and region spec

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




# References
* https://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.io.netcdf.netcdf_file.html

# Command line How to do
I think we'll need to a `__init__ == "__main__":` and run the command using
`python -m limited-area [options] grid.nc points.txt


## Packages we'll need
* Numpy
* Scipy-IO NetCDF - YES! It works!



# Points.txt syntax

Possible names for limited area
* limited_area
* area_def
* mpas_region
* regional_area_definition
* mpas_region_definition


# Possible Points Syntax

``` For all:
name: Region Name
type: circle/ellipse/square/custom
```

Circle with center and a radius
```
type: circle
center: lat, lon
radius: radius (with units: m/km/miles/feet)
```

Ellipse with left, right defined and then angle
```
type: ellipse
left: lat, lon
right: lat, lon
minor_axes_angle: angle (degrees)
```

Custom with a list a points. Point inside point needs to be spefieid with
`Point:`. Possible `inPoint`, `inside` notsure
```
type: custom
Point: Point inside point
lat1, lon1
lat2, lon2
lat3, lon3
lat4, lon4
lat5, lon5
```

Custom with pointing to a file. The file, or the points.txt file needs to
contain `point:` and at least two points(??)
```
type: custom
/path/to/file
```

Square with left, right, up down cords defined
```
type: square
left_lon: lon
right_lon: lon
up_lat: lat
down_lat: lat
```

Square with h and w and angle
```
type: square
h: height (units)
w: width (units)
angle: angle (default 0 - no ratation - rotate clockwise)
```

Optional arguments
```
radians - All latitude and longitude inputs are read in as radians
degrees - All lat and lon in degrees
meters/km - All distance units are read in as meters or km
borderSize: int - No less then 7 and produce warnings
algorithm: Follow-The-Line/Greedy/Dijkstra's
output: The output name of the mesh
```

## Current Points Syntax:

Circle: 
```
Name
circle
center_lat, center_lon, radius
```

Ellipse
```
Name
ellipse
lat_left, lon_left, lat_right, lon_right, semiminor_axes_angle
```

Custom
```
Name
Number of Points - Not Neccessary
lon1, lat1
lon2, lat2
...
lonN, latN
Point inside the region
```




