from __future__ import absolute_import, division, print_function
import argparse
import os
import sys

import numpy as np

from LimitedArea.mesh import MeshHandler, convert_lx
from LimitedArea.regionSpec import RegionSpec

class LimitedArea():
    ''' These could possible go into a settings.py file ?? '''
    BOUNDARY1 = 1
    BOUNDARY2 = 2
    BOUNDARY3 = 3
    BOUNDARY4 = 4
    BOUNDARY5 = 5
    BOUNDARY6 = 6
    BOUNDARY7 = 7
    INSIDE = 0
    UNMARKED = -1

    def __init__(self, 
                 mesh_files, # Mesh Filename
                 region,     # Region filename
                 regionFormat='points',  # Region Format
                 algorithm='follow',
                 *args, 
                 **kwargs):
        ''' Init function for Limited Area
    
        Check to see if mesh file exists and it is the correct type. Check to
        see that the region file exist and finally set the regionSpec to the
        requeste regionFormat


        Keyword arguments:
        mesh_files   -- Path to a valid MPAS Mesh file
        region       -- Path to a region specification - Can be either a points 
                        file or a shapeFile (Which is currently not 
                        implemented!)
        regionFormat -- The type of the region file. ie either 'points' or
                        'shapefile'

        Named Args TODO: Look up the format to do this
        DEBUG_VALUE -
        '''
        self.meshes = []

        self._DEBUG_VALUE = kwargs.get('DEBUG_VALUE', 0)
        self.algorithm = kwargs.get('algorithm', 'follow')

        # Check to see that all of the meshes exists and that they are netcdfs
        for mesh in mesh_files:
            if os.path.isfile(mesh):
                self.meshes.append(MeshHandler(mesh))

                if self._DEBUG_VALUE > 0: 
                    print("DEBUG: ", mesh, " is a valid NetCDF File\n")
            else:
                print("ERROR: Mesh file was not found", mesh)
                sys.exit(-1)

        # Check to see the points file exists and if it exists, then parse it
        # and see that is is specified correctly!
        if os.path.isfile(region):
            self.region_file = region
        else:
            print("ERROR: Region specification file was not found")
            sys.exit(-1)

        if regionFormat == 'points':
            self.regionSpec = RegionSpec(method=regionFormat)
            self.regionFormat = 'points'
        elif regionFormat == 'shape' .OR. regionFormat == 'shapeFile':
            self.regionSpec = RegionSpec(method=regionFormat)
            self.regionFormat = 'shape'
        else:
            raise NotImplementedError("REGION SPEC IS NOT IMPLMENTED "
                                      "- IMPEMTED IT!")

        if self.algorithm == 'follow':
            self.mark_boundry = self.follow_the_line
        elif self.algorithm == 'dijkstra':
            self.mark_boundry = self.dijkstra
        elif self.algorithm == 'greedy':
            self.mark_boundry = self.greedy 

    def gen_region(self, *args, **kwargs):
        ''' gen_region

        Flow Chart

        1. Call the region spec to generate `name, in_point, points`. Do this
           first. We don't want to open any of the mesh variables (which might be
           gigantic), so we can avoid unncesseary processing if the regionSpec
           encounters and error.
        
        '''

        # Call the regionSpec to generate `name, in_point, points`
        name, in_point, points = self.regionSpec.gen_spec(self.region_file)

        if self._DEBUG_VALUE > 0:
            print("DEBUG: Region Spec has been generated")
            print("DEBUG: Name: ", name)
            print("DEBUG: in_point: ", in_point)
            print("DEBUG: points: ", points)


        # For each mesh, mark the boundary
        
        for mesh in self.meshes:
            self.mark_boundry(mesh, in_point, points)
            self.flood_fill(mesh, in_point, boundary)
    

    def gen_output_filename(self, mesh_file, points):
        pass


    ''''''''''''''''''''''''''''''
    ''''''''''''''''''''''''''''''
    
    ''' Custom Algorithms '''

    def follow_the_line(self, mesh, in_point, points, *args, **kwargs):
        ''' Mark the nearest cell to each of the cords in points
        as a boundary cell.
        '''

        boundry_cells = []
        nCells = mesh.mesh.dimensions['nCells']
        latCell = mesh.mesh.variables['latCell'][:]
        lonCell = mesh.mesh.variables['lonCell'][:]
        sphere_radius = mesh.mesh.sphere_radius

        # Find the nearest cells to the list of given boundary points
        for i in range(0, len(points), 2):
            boundry_cells.append(mesh.nearest_cell(points[i],
                                                   points[i+1]))

        # Find the nearest cell to the inside point
        inside_point = mesh.nearest_cell(in_point[0], 
                                         in_point[1])

        print("DEBUG: Boundary Cells: ", boundry_cells)
        print("DEBUG: Inside point: ", inside_point)


        # Create the bdyMask fields
        # TODO: Update this with dtype,
        # TODO: Update this with order, C or F order in memory?
        bdyMaskCell = np.full(len(nCells), self.UNMARKED)
        bdyMaskEdge = np.full(len(nCells), self.UNMARKED)
        bdyMaskVertex = np.full(len(nCells), self.UNMARKED)
        

        # TODO: Test this with the enumerate Built-in function

        for i in range(len(boundry_cells)):
            source_cell = boundry_cells[i]
            target_cell = boundry_cells[i % len(boundry_cells)]

            xs, ys, zs = convert_lx(latCell[source_cell], 
                                    lonCell[source_cell], 
                                    sphere_radius)
            xt, yt, zt = convert_lx(latCell[target_cell],
                                    lonCell[target_cell],
                                    sphere_radius)
        
            print("DEBUG: Source Cell x, y, z: ", xs, ys, zs)
            print("DEBUG: Target Cell x, y, z: ", xt, yt, zt)

            cross = np.cross([xs, ys, zs], [xt, yt, zt])
            cross = cross / np.linalg.norm(cross) # Unit Vector

            
            
    def greedy(self, mesh, in_point, points, *args, **kwargs):
        pass

    def dijsktra(self, mesh, in_point, points, *args, **kwargs):
        pass
