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
$ # Running the script with no arguments will produce a usage output
$ create_region
Usage: create_region [-h] [-o OUTPUT] [-a ALGORITHM] [-v VERBOSE]
                     points grid [grid ...]
create_region: error: the following arguments are required: points, grid
```

**Note**: It may be necessary to install the dependencies for this
program if they have not been currently installed by you, or your administrator.
You can install all the dependencies for this repository by running  
`pip install -r requirements.txt`. This will install all the necessary 
dependencies needed to run this program.


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


## Notes on Creating Large Regions (nCells >= 2,000,000)<a name="large-Regions">

If the region you create is significantly large ( >= 2,000,000 grid cells) you
will need to change the NetCDF file version of the regional grid file. To do
this, you can change the `format` keyword argument of the `LimitedArea`
initialization call within the `create_region` script to `NETCDF3_64BIT_DATA`:

```
regional_area = LimitedArea(args.grid,
                            args.points,
                            format='NETCDF3_64BIT_DATA',
                            **kwargs)
```

## Notes on creating regions from .grid.nc files<a name="concave">

When creating region subsets from `grid.nc` files, it is important to take a
few things into account. Because of the way the init_atmosphere model 
interpolates some static data, it is possible to create incorrect meshes
when subsetting a global `grid.nc` file.

Creating a mesh that is either **concave** in overall shape, or spans
significant longitudes near the poles will cause init_atmosphere model
to incorrectly interpolate some static fields.

Thus, when subsetting a region of a `grid.nc` file, it is import to ensure
the subsetted region's shape is **convex** in shape and is significantly 
far enough away from the poles.

If you would like to have a convex region or region that spans many longitudes
it is possible to subvert this pitfal by creating regions from `static.nc` file.

# Points (pts) Syntax<a name="Points">

A number of pts syntax files are available in `docs/points-examples`.

The points syntax is a simple file format that can be used to specify a region.
In the future there will be a number of ways to specify different shapes with
the points syntax including, squares, circles, and ellipses, but currently the
only methods supported are the polygon and circle method.

Each method will differ slightly in the syntax used to describe a region, but a
number of keywords will be required by each method.

Each pts file will need to specify the following keywords followed by an
appropriate value:

 1. `Name:` - The desired name of your region. The specified keyword will be
              appended to the outputed region.
 2. `Type:` - The method for generating the region [custom/circle]
 3. `Point: `- A latitude, longitude coordinate separated by a comma of a point
               that is inside the desired region. **Note**: For the circle
               method, this point will be used as the center of the circle.


The value after `Name:` will be the name of the new regional mesh. 

If the `Name` within the pts file was set to be the following:
```
Name: region_name
```

Then the resulting regional MPAS grid will be named the following:
`region_name.10242.grid.nc` if the region was created from the
`x1.10242.grid.nc`.

## Polygon<a name="polygon">

A polygon points file would look like the following:
```
Name: Desired_name_of_region 
Type: custom        # Options are: [custom | circle | ellipse]
Point: lat, lon     # Point that inside the region
lat1, lon1          # List of points specifying the region
lat2, lon2
lat3, lon3
...
latN, LonN
```

For the polygon method, the value of `Type:` must be specified as `custom`.
The value of `Point:` must be a latitude, longitude point within the desired
region in degrees.

After the `Point` specification, any number of coordinate points (in degrees)
that define the desired boundary can be listed. Points should be listed in
counter-clockwise order.

An example polygon method for defining CONUS region (`conus.custom.pts`):
```
Name: conus
Type: custom
Point: 40.0, -100.0
50.0, -129.0
50.0, -65.0
20.0, -65.0
20.0, -129.0
```

**NOTE:** When creating regional meshes with the polygon method it is
neccessary to **not create concave** regions. Doing so will result in MPAS
init_atmosphere producing erronus results while interpolating static fields.
Please see the section on [Concave Regions](concave).

## Circle<a name="circle">

The circle method produces a regional circle subset of a global grid given a
center point (`Point`) and a radius from that point (`Radius`).

An example circle points specification would look like the following:
```
Name: my_circle
Type: circle
Point: -40.0, 105.5
Radius: 3000
```

Here `Point` is used to specify the center of the circle and `Radius` is used
to specify the radius of the circle (in kilometers). **NOTE**: The radius must
be larger then at least the smallest grid cell, else unexpected behavior will
occur.

An example circle method for defining a circle around Boulder, Colorado is:
```
Name: boulder
Type: circle
Point: 40.0, -105.5
radius: 4000
```

## Equatorial Channel

The equatorial channel method can be used to specify a channel parallel to the
equator, across longitudes and between two latitude points.

For instance, one can specify a region between the tropics as:

```
Name: tropics
Type: channel
Upper-lat: 23.43676
lower-lat: -23.43676
```

### **Caution on creating channels of `grid.nc` files**

Depending on the location of a specified channel, erroneous results can be
caused during some interpolation static fields by the init_atmosphere core.
This will cause erroneous results within result `static.nc` files and will
cause resulting files to be unusable.

Thus, it is **highly** recommended to only use the channel method upon
`static.nc` files. Doing so will avoid this problem and will also save time
interpolating the static fields on different regions.

# Reporting Bugs

If you encounter a bug and wish to report it, please do so on this Github repository's Issues page! Thank you!
