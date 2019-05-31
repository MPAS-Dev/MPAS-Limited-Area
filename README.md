# MPAS Limited-Area 

MPAS Limited-Area is a python tool that takes an MPAS global grid and produces
a regional area grid given a region specifications. Specifications can be
specified in a number of different ways, making limited-area extensible and
flexible.


# Download and Installing<a name="Installing"/>

To download this command line script, clone this repository into the location
of your choosing. To install it, add the base directory to your path i.e.:

```
$ git clone git@github.com:MiCurry/MPAS-Limited-Area.git
$ setenv PATH ${PATH}:/path/to/MPAS-Limited-Area
```

Then, run the `create_region` command-line script:
```Bash
$ # Running the script with no arguments will produce help output
$ create_region
Usage: create_region [-h] [-o OUTPUT] [-a ALGORITHM] [-v VERBOSE]
                     points grid [grid ...]
create_region: error: the following arguments are required: points, grid
```


# Running<a name="Running"/>

At any point in time, you can pass the `-h` or `--help` flag to `limited-area`
to generate a help message, usage statement and a list of options. The command
line usage is:

``` bash
$ create_region [options] points grid
```

Where `points` is a points specification file, and where `grid` is a global
MPAS NetCDF grid file. You'll find example points file in the points-example
directory within the docs directory of this repository and below in the *Points
Syntax* section.

# Points Syntax<a name="Points">


The points syntax is a simple file format that can be used to specify a region.
In the future there will be a number of ways to specify different shapes with
the points syntax including, squares, circles, and ellipses, but currently the
method is the polygon method.

A polygon points file would look like the following:
```
Name: Desired_name_of_region 
Type: custom    
Point: lat, lon 
lat1, lon1
lat2, lon2
lat3, lon3
...
latN, LonN
```

Where the value after `Name:` will be the name of the new regional mesh. The
regional filename will follow the convention: `region_name.10242.grid.nc` if
the region was created from the `x1.10242.grid.nc`.

For the polygon method, the value of `Type:` must be specified as `custom`.
The value of `Point:` must be a latitude, longitude point within the desired
region in degrees.

After the `Point` specification, any number of coordinate points (in degrees)
that define the desired boundary can be listed. Points should be listed in
counter-clockwise order.

## CONUS Points Example

Define a polygon method by defining a new points file ( for instance,
`conus.custom.pts` ) as the following:

```
Name: conus
Type: custom
Point: 40.0, -100.0
50.0, -129.0
50.0, -65.0
20.0, -65.0
20.0, -129.0
```

