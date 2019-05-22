# MPAS Limited-Area 

MPAS Limited-Area is a python tool that takes an MPAS global grid and produces
a regional area grid given a region specifications. Specifications can be
specified in a number of different ways, making limited-area extensible and
flexible.

# Table of Contents

1. [Downloading & Installing](#Installing)
2.  [Running](#Running)

# Download & Installing<a name="Installing"/>

To download and install this command line script, we will clone this
repository, and then we will use the `setupt.py` file to install the script to
a location of our choosing.

0. Insure that we have the correct dependencies:

Install the python dependencies, if they are not already installed. If pip
responds with a 'permission denied', you may need to use the '--user' to
install these.
```
> pip install numpy
> pip install netcdf4
```
**NOTE:** For installing `netcdf4`, you'll need to insure you have the correct
binary requirements (HDF5 C libraries installed, netCDF-4 C Library install,
etc). Please see <http://unidata.github.io/netcdf4-python/netCDF4/index.html>
for a complete list of dependencies for Python NetCDF4 as well as install
instructions.

1. Clone this repository
```
> git clone https://github.com/MiCurry/MPAS-Limited-Area.git
```

2. Make the file `limited-area` runnable: 
```
> chmod 750 limited-area
```

Now add this directory to your `PATH` environment variable

```
# csh shells
> setenv PATH ${PATH}:/path/to/MPAS-Limited-Area
```

```
# Bash shells
> export PATH=${PATH}:/path/to/MPAS-Limited-Area
```

3. Test the `limited-area` is on your path and is runnable by running the
`limited-area`:

```
> limited-area
```

If it runs successfully you should get a small help message:
```
usage: limited-area [-h] [-o OUTPUT] [-a ALGORITHM] [-v VERBOSE]
                    points grid [grid ...]
                    limited-area: error: the following arguments are required:
                    points, grid
```

# Running<a name="Running"/>

At any point in time, you can pass the `-h` or `--help` flag to `limited-area`
to generate a help message and usage statement. The command line usage is:

```
limited-area [options] points grid
```

Where `points` is a points specification file, and where `grid` is a global
MPAS NetCDF grid file. You'll find example points file in the points-example
directory within the docs directory of this repository.
