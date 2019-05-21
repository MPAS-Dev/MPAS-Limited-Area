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

If we have sudo privileges on your machine, you can install `scipy`, `numpy`
and `netcdf4` using pip as one normally would:
```
pip install numpy
pip install scipy
pip install netcdf4
```
If you don't have sudo privileges, then you can install the same packages with
the `--user` option. i.e. `pip install numpy --user` etc.

**NOTE:** For installing `netcdf4`, you'll need to insure you have the correct
binary requirements (HDF5 C libraries installed, netCDF-4 C Library install,
etc). Please see <http://unidata.github.io/netcdf4-python/netCDF4/index.html>
for a complete list of dependencies for Python NetCDF4 as well as install
instructions.


1. Git clone this repository
```
git clone git@github.com:MiCurry/MPAS-Limited-Area.git
```

2. Install the limited-area script to a location of your choosing using the
following command and specifying the location with the `--home=<dir>` option.
```
python setup.pt install --home=/path/to/install/dir
```

3. If the location of the installation is not in your path, then add it to your
path:
```
export PATH=${PATH}:/path/to/install/dir
```

4. Test to see if the installation ran correctly by running:
```
limited-area
```

# Running<a name="Running"/>

At any point in time, you can pass the `-h` or `--help` flag to `limited-area`
to generate a help message and usage statement. The command line usage is:

```
limited-area [options] points grid
```

Where `points` is a points specification file, and where `grid` is a global
MPAS NetCDF grid file.

