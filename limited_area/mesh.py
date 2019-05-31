from __future__ import absolute_import, division, print_function
import sys
import os

import numpy as np
from netCDF4 import Dataset

""" mesh.py - Handle NetCDF file operations as well as calculations
upon on MPAS grid."""

class MeshHandler:
    """ Handle the operations related to NetCDF/MPAS grids. """

    def __init__(self, fname, mode, *args, **kwargs):
        """ """
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
        """ Create and open a new NetCDF file with name fname and access mode mode """
        try:
            self.mesh = Dataset(fname, mode)
            return
        except:
            print("ERROR: There was a problem creating the file ", fname)
            sys.exit(-1)


    def check_file(self, fname):
        """ Check to see that fname exists and it is a valid NetCDF file """
        if os.path.isfile(fname):
            try:
                self.mesh = Dataset(fname, 'r')
                if self._DEBUG_ > 2:
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
        """ Print all the dimension names for this mesh"""
        print(self.mesh.dimensions.keys())


    def print_all_variables(self):
        """ Print all the variable names for this mesh"""
        print(self.mesh.variables.keys())


    def nearest_cell(self, lat, lon):
        """ Find the nearest cell of this mesh to lat and lon

        lat - Latitude - Radians
        lon - Longitude - Radians
        """
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

    def create_graph_file(self, graphFname):
        nCells = self.mesh.dimensions['nCells'].size
        nEdges = self.mesh.dimensions['nEdges'].size

        nEdgesOnCell = self.mesh.variables['nEdgesOnCell'][:]
        cellsOnCell = self.mesh.variables['cellsOnCell'][:]
        cellsOnEdge = self.mesh.variables['cellsOnEdge'][:]
        
        nEdgesInterior = 0
        for i in range(nEdges):
            if cellsOnEdge[i,0] > 0 and cellsOnEdge[i,1] > 0:
                nEdgesInterior = nEdgesInterior + 1


        with open(graphFname, 'w') as f:
            f.write(repr(nCells)+' '+repr(nEdgesInterior)+'\n')
            for i in range(nCells):
                for j in range(nEdgesOnCell[i]):
                    if (cellsOnCell[i,j] > 0):
                        f.write(repr(cellsOnCell[i,j])+' ')
                f.write('\n')
            
        print(graphFname)

    def subset_fields(self, 
                      regionalFname, 
                      bdyMaskCell,
                      bdyMaskEdge,
                      bdyMaskVertex,
                      inside,
                      unmarked,
                      *args, 
                      **kwargs):
        """ Subset the current mesh and return a new regional mesh with
        subsetted fields 
        
        regionalFname -- Desired filename for the regional subset
        bdyMaskCell   -- Global mesh mask denoting regional cells
        bdyMaskEdge   -- Global mesh mask denoting regional edges
        bdyMaskVertex -- Global mesh mask denoting regional vertices
        inside        -- The integer value that was used to mark the 
                         cells, edges, vertices as being 'inside' the 
                         regional within the bdyMasks
        unmarked      -- The integer value that was used to mark cells,
                         edges, vertices as being 'outside' of the regional
                         mesh.
        """

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

        # Check to see the user didn't mess specifying the region. If 
        # len(bdyIndexToCellIDs) == nCells, then the specification was probably not
        # specified correctly
        if len(bdyIndexToCellIDs) == nCells:
            print("ERROR: The number of Cells in the specified region ",
                  "(", len(bdyIndexToCellIDs), ")")
            print("ERROR: appears to be equal number of cells in the global mesh",
                  "(", nCells, ")")
            print("ERROR: which means there was perhaps a problem in specifying the")
            print("ERROR: region. Please insure your region specification is correct")
            sys.exit(-1)
        

        # Create a new grid
        region = MeshHandler(regionalFname, 'w', *args, **kwargs)

        # Dimensions - Create dimensions
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
        
        # Variables - Create Variables
        for var in self.mesh.variables:
            region.mesh.createVariable(var, self.mesh.variables[var].dtype,
                                            self.mesh.variables[var].dimensions)

        # Subset global variables into the regional mesh and write them
        # to the regional mesh - reindexing if neccessary
        for var in self.mesh.variables:
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


        # Make boundary Mask's between 0 and the number of speficified relaxtion
        # layers
        region.mesh.variables['bdyMaskCell'][:] = bdyMaskCell[bdyMaskCell != 0] - 1 
        region.mesh.variables['bdyMaskEdge'][:] = bdyMaskEdge[bdyMaskEdge != 0] - 1
        region.mesh.variables['bdyMaskVertex'][:] = bdyMaskVertex[bdyMaskVertex != 0] - 1

        return region

    def copy_global_attributes(self, region):
        """ Copy the global attributes into the regional mesh, but not 'np' """
        region.mesh.np = region.mesh.dimensions['nCells'].size

        region.mesh.on_a_sphere = self.mesh.on_a_sphere
        region.mesh.sphere_radius = self.mesh.sphere_radius
        region.mesh.n_scvt_iterations = self.mesh.n_scvt_iterations
        region.mesh.eps = self.mesh.eps
        region.mesh.Convergence = self.mesh.Convergence




def reindex_field(field, mmap):
    """ If field[i] is in mmap, then reindex it with the index of 
    mmap where it equals field[i] """
    print(' ... Reindexing Field ... ', field.shape, mmap.shape, end=' ... ', flush=True)
    field = field.flatten()
    for i in range(len(field)):
        field[i] = binary_search(mmap, field[i])
    
    print(' Done!')
    return field


def binary_search(arr, x):
    """ Search for x in arr and return the location of x in arr or
    return 0 if it is not found. """
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


def latlon_to_xyz(lat, lon, radius):
    """ Calculate and return x, y, z cordinations of lat, lon on the sphere that has
    radius, radius.
    lat - Latitude
    lon - Longitutde
    radius - Radius of sphere
    """
    z = radius * np.sin(lat)
    x = radius * np.cos(lon)
    y = radius * np.sin(lon) * np.cos(lat)

    return np.array([x, y, z])


def sphere_distance(lat1, lon1, lat2, lon2, radius, **kwargs):
    """ Calculate the sphere distance between point1 and point2. 

    lat1 - Float - Radians - -pi:pi
    lon1 - Float - Radians - 0:2*pi
    lat2 - Float - Radians - -pi:pi
    lon2 - Float - Radians - 0:2*pi
    radius - Radius of the earth (or sphere) - Units can be ignored

    """ 
    return (2 * radius * np.arcsin(
                         np.sqrt(
                         np.sin(0.5 * (lat2 - lat1))**2
                       + np.cos(lat1) 
                       * np.cos(lat2) 
                       * np.sin(0.5 * (lon2 - lon1))**2)))
