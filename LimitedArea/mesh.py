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
        nCells = self.mesh.dimensions['nCells']
        latCells = np.array(self.mesh.variables['latCell'])
        lonCells = np.array(self.mesh.variables['lonCell'])
        nEdgesOnCell = np.array(self.mesh.variables['nEdgesOnCell'])
        cellsOnCell = np.array(self.mesh.variables['cellsOnCell'])
        sphere_radius = self.mesh.sphere_radius

        nearest_cell = 1 # Start at the first cell
        current_cell = -1


        print("DEBUG: Nearest_cell Lat: ", lat, " Lon: ", lon)

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
                iCell = cellsOnCell[current_cell, edges]
                if (iCell <= nCells.size):
                    iDistance = sphere_distance(latCells[iCell],
                                                lonCells[iCell],
                                                lat,
                                                lon,
                                                sphere_radius)

                    if (iDistance < nearest_distance):
                        nearest_cell = iCell
                        nearest_distance = iDistance
    
        if self._DEBUG_ > 5:
            print("DEBUG: nearest_cell latLon: ", nearest_cell,
                                                  latCells[nearest_cell],
                                                  lonCells[nearest_cell])
            print("DEBUG: Nearest Cell: ", nearest_cell)

        return nearest_cell


    def subset_fields(self, 
                      regionalFname, 
                      bdyMaskCell,
                      bdyMaskEdge,
                      bdyMaskVertex,
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


       # print("bdyMaskCell: ", bdyMaskCell, len(bdyMaskCell))
       # input('-')
       # print("bdyMaskEdge: ", bdyMaskVertex, len(bdyMaskEdge))
       # input('-')
       # print("bdyMaskVertex: ", bdyMaskVertex, len(bdyMaskVertex))
       # input('-')

        bdyIndexToCellIDs = indexToCellIDs[np.where(bdyMaskCell != 0)]
        bdyIndexToEdgeIDs = indexToEdgeIDs[np.where(bdyMaskEdge != 0)]
        bdyIndexToVertexIDs = indexToVertexIDs[np.where(bdyMaskVertex != 0)]

       # print("bdyIndexToCellIDs: ", bdyIndexToCellIDs, " size: ", len(bdyIndexToCellIDs))
       # input("-")
       # print("bdyIndexToEdgeIDs: ", bdyIndexToEdgeIDs, " size: ", len(bdyIndexToEdgeIDs))
       # input("-")
       # print("bdyIndexToVertexIDs: ", bdyIndexToVertexIDs, " size: ", len(bdyIndexToVertexIDs))
       # input("-")


        # Create a new grid
        region = MeshHandler(regionalFname, 'w', *args, **kwargs)

        ''' Dimensions '''
        for dim in self.mesh.dimensions:
            if dim == 'nCells':
                region.mesh.createDimension(dim, len(bdyIndexToCellIDs))
            elif dim == 'nEdges':
                region.mesh.createDimension(dim, len(bdyIndexToEdgeIDs))
            elif dim == 'nVertices':
                region.mesh.createDimension(dim, len(bdyIndexToVertexIDs))
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
                                     reindex_field(self.mesh.variables[var][bdyIndexToCellIDs],
                                                              bdyIndexToEdgeIDs)
            elif var == 'verticesOnCell':
                region.mesh.variables[var][:] = \
                                    reindex_field(self.mesh.variables[var][bdyIndexToCellIDs],
                                                  bdyIndexToVertexIDs)
            elif var == 'cellsOnCell':
                region.mesh.variables[var][:] = \
                                    reindex_field(self.mesh.variables[var][bdyIndexToCellIDs],
                                                  bdyIndexToCellIDs)
            # Vertices
            elif var == 'edgesOnVertex':
                region.mesh.variables[var][:] = \
                                    reindex_field(self.mesh.variables[var][bdyIndexToVertexIDs],
                                                  bdyIndexToEdgeIDs)
            elif var == 'cellsOnVertex':
                region.mesh.variables[var][:] = \
                                    reindex_field(self.mesh.variables[var][bdyIndexToVertexIDs],
                                                  bdyIndexToCellIDs)
            # Edges
            elif var == 'verticesOnEdge':
                region.mesh.variables[var][:] = \
                                    reindex_field(self.mesh.variables[var][bdyIndexToEdgeIDs],
                                                  bdyIndexToVertexIDs)
            elif var == 'cellsOnEdge':
                region.mesh.variables[var][:] = \
                                    reindex_field(self.mesh.variables[var][bdyIndexToEdgeIDs],
                                                  bdyIndexToCellIDs)
            elif var == 'edgesOnEdge':
                region.mesh.variables[var][:] = \
                                    reindex_field(self.mesh.variables[var][bdyIndexToEdgeIDs],
                                                  bdyIndexToEdgeIDs)
                

            elif 'nCells' in self.mesh.variables[var].dimensions:
                print(" ")
                if var in indexingFields:
                    region.mesh.variables[var][:] = \
                                        reindex_field(self.mesh.variables[var][bdyIndexToCellIDs],
                                                      bdyIndexToCellIDs)
                else:
                    region.mesh.variables[var][:] = self.mesh.variables[var][bdyIndexToCellIDs]
            elif 'nEdges' in self.mesh.variables[var].dimensions:
                print(" ")
                if var in indexingFields:
                    region.mesh.variables[var][:] = \
                                        reindex_field(self.mesh.variables[var][bdyIndexToEdgeIDs],
                                                      bdyIndexToEdgeIDs)
                else:
                    region.mesh.variables[var][:] = self.mesh.variables[var][bdyIndexToEdgeIDs]
            elif 'nVertices' in self.mesh.variables[var].dimensions:
                print(" ")
                if var in indexingFields:
                    region.mesh.variables[var][:] = \
                                        reindex_field(self.mesh.variables[var][bdyIndexToVertexIDs],
                                                      bdyIndexToVertexIDs)
                else:
                    region.mesh.variables[var][:] = self.mesh.variables[var][bdyIndexToVertexIDs]

        return region


def reindex_field(field, index):
    newField = np.searchsorted(index, field, side='right')
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
                         np.sin(0.5 * (lat2 - lat1)**2 
                       + np.cos(lat1) 
                       * np.cos(lat2) 
                       * np.sin(0.5 * (lon2 - lon1))**2))))
