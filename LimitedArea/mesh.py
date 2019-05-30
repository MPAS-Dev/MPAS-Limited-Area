''' NetCDF Handler '''
from __future__ import absolute_import, division, print_function
import sys
import os

import numpy as np
from netCDF4 import Dataset


class MeshHandler:
    def __init__(self, fname, mode, *args, **kwargs):
        self._DEBUG_ = kwargs.get('DEBUG', 0)
        self.fname = fname

        if mode == 'r':
            if self.check_file(fname):
                return
            else:
                sys.exit(-1)
        elif mode == 'w':
            self.create_file(fname, mode)

    def create_file(self, fname, mode):
        try:
            self.mesh = Dataset(fname, mode)
            return
        except:
            print("ERROR: There was a problem creating the file ", fname)
            sys.exit(-1)

    def check_file(self, fname):
        if os.path.isfile(fname):
            try:
                self.mesh = Dataset(fname, 'r')
                if self._DEBUG_ > 1:
                    print("DEBUG: Mesh's dimensions: ", fname)
                    self.print_all_dimensions()
                    print("DEBUG: Mesh's variables: ", fname)
                    self.print_all_variables()
                return True
            except OSError as E: 
                print("ERROR: ", E)
                print("ERROR: This file was not a valid NetCDF file")
                sys.exit(-1)
        else:
            print("ERROR: This file did not exist!")
            return False


    def print_all_dimensions(self):
        print(self.mesh.dimensions.keys())

    def print_all_variables(self):
        print(self.mesh.variables.keys())

    def nearest_cell(self, lat, lon):
        ''' nearest_cell

        Find the nearest cell to lat and lon

        lat - Latitude - In degrees - -90:90
        lon - Longitude  - In degrees - -180:180

        '''
        # We will start at cell 0
        nCells = self.mesh.dimensions['nCells'].size
        latCells = np.array(self.mesh.variables['latCell'])
        lonCells = np.array(self.mesh.variables['lonCell'])
        nEdgesOnCell = np.array(self.mesh.variables['nEdgesOnCell'])
        cellsOnCell = np.array(self.mesh.variables['cellsOnCell'])
        sphere_radius = self.mesh.sphere_radius

        nearest_cell = 0 # Start at the first cell
        current_cell = -1

        while (nearest_cell != current_cell):
            current_cell = nearest_cell
            current_distance = sphere_distance(latCells[current_cell],
                                               lonCells[current_cell],
                                               lat,
                                               lon,
                                               sphere_radius)
            
            nearest_cell = current_cell
            nearest_distance = current_distance
            
            for edges in range(nEdgesOnCell[current_cell]):
                iCell = cellsOnCell[current_cell, edges] - 1
                if (iCell <= nCells):
                    iDistance = sphere_distance(latCells[iCell],
                                                lonCells[iCell],
                                                lat,
                                                lon,
                                                sphere_radius)

                    if (iDistance <= nearest_distance):
                        nearest_cell = iCell
                        nearest_distance = iDistance

    
        if self._DEBUG_ > 3:
            print("DEBUG: nearest_cell latLon: ", nearest_cell, '\t',
                                                  latCells[nearest_cell],
                                                  lonCells[nearest_cell],
                  ' Given lat lon: ', lat, lon)


        return nearest_cell


    def subset_fields(self, 
                      regionalFname, 
                      bdyMaskCell,
                      bdyMaskEdge,
                      bdyMaskVertex,
                      inside,
                      unmarked,
                      *args, 
                      **kwargs):
        ''' From the current mesh, subset all the fields into a new mesh,
        regionalFname  '''

        # Don't pass on DEBUG to the regional mess - tone down output
        kwargs.pop('DEBUG')

        indexingFields = [ 'indexToCellID', 'indexToEdgeID', 'indexToVertexID' ]
                            
        nCells = self.mesh.dimensions['nCells'].size
        indexToCellIDs = self.mesh.variables['indexToCellID']
        indexToEdgeIDs = self.mesh.variables['indexToEdgeID']
        indexToVertexIDs = self.mesh.variables['indexToVertexID']

        bdyIndexToCellIDs = indexToCellIDs[np.where(bdyMaskCell != unmarked)]
        bdyIndexToEdgeIDs = indexToEdgeIDs[np.where(bdyMaskEdge != unmarked)]
        bdyIndexToVertexIDs = indexToVertexIDs[np.where(bdyMaskVertex != unmarked)]


        # Create a new grid
        region = MeshHandler(regionalFname, 'w', *args, **kwargs)

        ''' Dimensions '''
        for dim in self.mesh.dimensions:
            if dim == 'nCells':
                region.mesh.createDimension(dim, 
                                            len(bdyIndexToCellIDs))
            elif dim == 'nEdges':
                region.mesh.createDimension(dim, 
                                            len(bdyIndexToEdgeIDs))
            elif dim == 'nVertices':
                region.mesh.createDimension(dim, 
                                            len(bdyIndexToVertexIDs))
            else:
                region.mesh.createDimension(dim,
                                            self.mesh.dimensions[dim].size)
        
        ''' Variables '''
        # Create variables
        for var in self.mesh.variables:
            region.mesh.createVariable(var, self.mesh.variables[var].dtype,
                                            self.mesh.variables[var].dimensions)

        for var in self.mesh.variables: # Subset and write variables
            print("Copying variable: ", var, end=' ', flush=True)

            # Cells
            if var == 'edgesOnCell':
                region.mesh.variables[var][:] = \
                                     reindex_field(self.mesh.variables[var][bdyIndexToCellIDs-1],
                                                              bdyIndexToEdgeIDs)
            elif var == 'verticesOnCell':
                region.mesh.variables[var][:] = \
                                    reindex_field(self.mesh.variables[var][bdyIndexToCellIDs-1],
                                                  bdyIndexToVertexIDs)
            elif var == 'cellsOnCell':
                region.mesh.variables[var][:] = \
                                    reindex_field(self.mesh.variables[var][bdyIndexToCellIDs-1],
                                                  bdyIndexToCellIDs)
            # Vertices
            elif var == 'edgesOnVertex':
                region.mesh.variables[var][:] = \
                                    reindex_field(self.mesh.variables[var][bdyIndexToVertexIDs-1],
                                                  bdyIndexToEdgeIDs)
            elif var == 'cellsOnVertex':
                region.mesh.variables[var][:] = \
                                    reindex_field(self.mesh.variables[var][bdyIndexToVertexIDs-1],
                                                  bdyIndexToCellIDs)
            # Edges
            elif var == 'verticesOnEdge':
                region.mesh.variables[var][:] = \
                                    reindex_field(self.mesh.variables[var][bdyIndexToEdgeIDs-1],
                                                  bdyIndexToVertexIDs)
            elif var == 'cellsOnEdge':
                region.mesh.variables[var][:] = \
                                    reindex_field(self.mesh.variables[var][bdyIndexToEdgeIDs-1],
                                                  bdyIndexToCellIDs)
            elif var == 'edgesOnEdge':
                region.mesh.variables[var][:] = \
                                    reindex_field(self.mesh.variables[var][bdyIndexToEdgeIDs-1],
                                                  bdyIndexToEdgeIDs)
                
            elif var == 'indexToCellID':
                print(" ")
                region.mesh.variables[var][:] = np.arange(1, len(bdyIndexToCellIDs) + 1)
            elif var == 'indexToEdgeID':
                print(" ")
                region.mesh.variables[var][:] = np.arange(1, len(bdyIndexToEdgeIDs) + 1)
            elif var == 'indexToVertexID':
                print(" ")
                region.mesh.variables[var][:] = np.arange(1, len(bdyIndexToVertexIDs) + 1)
            elif 'nCells' in self.mesh.variables[var].dimensions:
                print(" ")
                if var in indexingFields:
                    region.mesh.variables[var][:] = \
                                        reindex_field(self.mesh.variables[var][bdyIndexToCellIDs-1],
                                                      bdyIndexToCellIDs)
                else:
                    region.mesh.variables[var][:] = self.mesh.variables[var][bdyIndexToCellIDs-1]
            elif 'nEdges' in self.mesh.variables[var].dimensions:
                print(" ")
                if var in indexingFields:
                    region.mesh.variables[var][:] = \
                                        reindex_field(self.mesh.variables[var][bdyIndexToEdgeIDs-1],
                                                      bdyIndexToEdgeIDs)
                else:
                    region.mesh.variables[var][:] = self.mesh.variables[var][bdyIndexToEdgeIDs-1]
            elif 'nVertices' in self.mesh.variables[var].dimensions:
                print(" ")
                if var in indexingFields:
                    region.mesh.variables[var][:] = \
                                        reindex_field(self.mesh.variables[var][bdyIndexToVertexIDs-1],
                                                      bdyIndexToVertexIDs)
                else:
                    region.mesh.variables[var][:] = self.mesh.variables[var][bdyIndexToVertexIDs-1]


        region.mesh.createVariable('bdyMaskCell', 'i4', ('nCells',))
        region.mesh.createVariable('bdyMaskEdge', 'i4', ('nEdges',))
        region.mesh.createVariable('bdyMaskVertex', 'i4', ('nVertices')) 

        region.mesh.variables['bdyMaskCell'][:] = bdyMaskCell[bdyMaskCell != 0] - 1 
        region.mesh.variables['bdyMaskEdge'][:] = bdyMaskEdge[bdyMaskEdge != 0] - 1
        region.mesh.variables['bdyMaskVertex'][:] = bdyMaskVertex[bdyMaskVertex != 0] - 1

        return region


def reindex_field_slow(field, mmap):
    print(' ... Reindexing Field ... ', field.shape, mmap.shape, end=' ... ', flush=True)

    shape = field.shape
    field = field.ravel()
    newField = np.zeros(len(field))

    for f in range(len(field)):
        for m in range(len(mmap)):
            if field[f] == mmap[m]: 
                newField[f] = m + 1
                break

    print(' Done!')
    return newField


def reindex_field(field, mmap):
    print(' ... Reindexing Field ... ', field.shape, mmap.shape, end=' ... ', flush=True)
    field = field.flatten()
    for i in range(len(field)):
        field[i] = binary_search(mmap, field[i])
    
    print(' Done!')
    return field

def binary_search(arr, x):
    l = 0
    u = len(arr) - 1
    k = (l + u) // 2

    while u >= l:
        if arr[k] == x:
            return k + 1
        elif arr[k] < x:
            l = k + 1
            k = (l + u) // 2
        else:
            u = k - 1
            k = (l + u) // 2

    return 0


def _reindex_field(field, mmap):
    print(' ... Reindexing Field ... ', field.shape, mmap.shape)

    shape = field.shape
    field = field.ravel()

    np.set_printoptions(threshold=np.inf)
    print('mmap: ', mmap)

    for i in range(len(field)):
        numCuts = 0
        cut = int(np.floor(mmap.shape[0]/2))

        while True:
            numCuts += 1
            if numCuts >= int(mmap.shape[0]): 
                print('Did not find: ', field[i])
                field[i] = 0
                input('-')
                break
            elif cut > mmap.shape[0]:
                print('Did not find: ', field[i])
                field[i] = 0
                input('-')
                break
            elif field[i] == mmap[cut]:
                field[i] = cut + 1
                break
            elif field[i] > mmap[cut]:
                cut += int(np.floor(cut/2))
            elif field[i] < mmap[cut]:
                cut = int(np.floor(cut/2))


    field = np.reshape(field, shape)

    return field


def reindex_field_basic(field, mmap):
    print('... Reindexing a basic field ...', field.shape, mmap.shape)
    newField = np.searchsorted(field, mmap) + 1
    return newField


def copy_mesh(inMesh):
    ''' Copy `inMesh` and then return a new mesh class that points to the copy
    '''
    pass


def latlon_to_xyz(lat, lon, radius):
    ''' Calculate the x, y, z cordinations of lat, lon on the sphere that has
    radius, radius.
    lat -
    lon -
    radius -
    TODO: Update this doc
    '''
    z = radius * np.sin(lat)
    x = radius * np.cos(lon)
    y = radius * np.sin(lon) * np.cos(lat)

    return np.array([x, y, z])


def sphere_distance(lat1, lon1, lat2, lon2, radius, **kwargs):
    ''' Calculate the sphere distance between point1 and point2. 

    lat1 - Float - Radians - -pi:pi
    lon1 - Float - Radians - 0:2*pi
    lat2 - Float - Radians - -pi:pi
    lon2 - Float - Radians - 0:2*pi
    radius - Radius of the earth (or sphere) - Units can be ignored

    TODO: Update doc with return value
    '''


    return (2 * radius * np.arcsin(
                         np.sqrt(
                         np.sin(0.5 * (lat2 - lat1))**2
                       + np.cos(lat1) 
                       * np.cos(lat2) 
                       * np.sin(0.5 * (lon2 - lon1))**2)))
