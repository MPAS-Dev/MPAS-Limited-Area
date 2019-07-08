from limited_area.mesh import MeshHandler
import sys
''' For each cell, find the lat and lon, use mesh.nearest_neighbor to find the
nn, and then brute force the nearest neighbor. '''
def nn_test(meshFile):
    mesh = MeshHandler(meshFile, 'r')
    
    nCells = mesh.mesh.dimensions['nCells'].size
    latCell = mesh.mesh.variables['latCell']
    lonCell = mesh.mesh.variables['lonCell']

    cell = 0
    for i in range(nCells):
        nn = mesh.nearest_cell(latCell[i], lonCell[i])
       
        if nn != i:
            print("FAILED!", 
                  "Given coord: ", latCell[i], lonCell[i], i,
                  "nc cord:", latCell[nn], lonCell[nn], nn)
            sys.exit(-1)
        else:
            print("PASSED!", 
                  "Given coord: ", latCell[i], lonCell[i],
                  "nc cord:", latCell[nn], lonCell[nn], nn)
