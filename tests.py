import argparse
import os

from tests.nearest_cell_test import nn_test
from tests.test_meshes import test_create_region


""" Test a number of meshes to ensure that they work correctly """
grids = ['x1.2562.grid.nc', 
         #'x1.4002.grid.nc',
          'x1.10242.grid.nc',
          'x1.40962.grid.nc',
          'x1.163842.grid.nc']

statics = ['x1.2562.static.nc', 
           'x1.10242.static.nc',
           'x1.40962.static.nc',
           'x1.655362.static.nc']

regions = ['conus.custom.pts',
           'india.circle.pts',
           'japan.ellipse.pts',
           'pnw.custom.pts',
           'tropics.channel.pts']

def test_grids():
    gridDir = '/glade/work/mcurry/meshes2/'
    regionDir = '/glade/u/home/mcurry/mpas/MPAS-Limited-Area/docs/points-examples'
    
    print("-------------------------------")
    print("Testing all example regions on grid files:")
    print("")

    for grid in grids:
        for region in regions:
            test_create_region(gridDir, 
                               grid, 
                               os.path.join(regionDir, region))
            print("")

def test_statics():
    staticDir = '/glade/work/mcurry/static_runs/'
    regionDir = '/glade/u/home/mcurry/mpas/MPAS-Limited-Area/docs/points-examples'

    print("-------------------------------")
    print("Testing all example regions on static files:")
    print("")

    for static in statics:
        for region in regions:
            test_create_region(staticDir, 
                               static, 
                               os.path.join(regionDir, region))
            print("")

    print("Tests complete!")
    print("")


if __name__ == "__main__":
    gridDir = '/glade/work/mcurry/meshes/'
    staticDir = '/glade/work/mcurry/static_runs/'
    regionDir = ''

    description = ( "Test launcher for MPAS-Limited-Area" )


    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('test',
                        help='Test name or suite name',
                        type=str)
    parser.add_argument('-r', '--region',
                        help='Test name or suite name',
                        type=str)

    args = parser.parse_args()


    print(args)

    if args.test == 'nn_test':
        meshFile= '/glade/work/mcurry/meshes/x1.2562.grid.nc'
        nn_test(meshFile)
    if args.test == 'grid':
        test_grids()
    if args.test == 'static':
        test_statics()
    if args.test == 'all':
        test_grids()
        test_statics()
