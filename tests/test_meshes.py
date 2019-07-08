import sys
import os

from limited_area.limited_area import LimitedArea
from netCDF4 import Dataset

import traceback

def test_create_region(meshDir, mesh, regionSpec, DEBUG=0):

    print( 'TEST_CREATE_REGION' )
    print( 'TEST_CREATE_REGION: mesh: ', os.path.join(meshDir, mesh))
    print( 'TEST_CREATE_REGION: region: ', regionSpec)
    print("")

    regional_area = LimitedArea(os.path.join(meshDir, mesh),
                               regionSpec,
                               DEBUG=DEBUG)

    try:
        region, regionGraph = regional_area.gen_region(DEBUG=DEBUG)
    except Exception as e:
        print("ERROR:", str(e)) 
        traceback.print_exc()
        print("\nFAILED ----------------------------------------------------------------------")
        print("TEST_CREATE_REGION: Test FAILED for grid ", mesh, " on region ", regionSpec)
        print("FAILED ----------------------------------------------------------------------\n")
        sys.exit(-1)

    if not os.path.isfile(region):
        print("\nFAILED ----------------------------------------------------------------------")
        print("TEST_CREATE_REGION: ERROR - region file: ", region," does not exist!")
        print("FAILED ----------------------------------------------------------------------\n")
        sys.exit(-1)

    if not os.path.isfile(regionGraph):
        print("\nFAILED ----------------------------------------------------------------------")
        print("TEST_CREATE_REGION: ERROR - region graph file: ", regionGraph," does not exist!")
        print("FAILED ----------------------------------------------------------------------\n")
        sys.exit(-1)

    try:
        f = Dataset(region)
        f.close()
    except:
        print(sys.exc_info()[0]) 
        print("\nFAILED ----------------------------------------------------------------------")
        print("TEST_CREATE_REGION: ERROR - region, ", region," could not be open by NetCDF4 Python!")
        print("FAILED ----------------------------------------------------------------------\n")

    del(regional_area) 

    print("")
    print("\nPASSED ----------------------------------------------------------------------")
    print("TEST_CREATE_REGION: All tests PASSED for ", mesh, " on region ", region)
    print("PASSED ----------------------------------------------------------------------\n")
    print("")
    return
