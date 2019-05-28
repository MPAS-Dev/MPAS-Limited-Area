from __future__ import absolute_import, division, print_function
import argparse
import os
import sys

import numpy as np

from LimitedArea.mesh import MeshHandler
from LimitedArea.mesh import latlon_to_xyz
from LimitedArea.mesh import sphere_distance
from LimitedArea.regionSpec import RegionSpec

class LimitedArea():
    ''' These could possible go into a settings.py file ?? '''
    num_boundary_layers = 7
    INSIDE = 1
    UNMARKED = 0

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

        ''' Named arguments '''
        self._DEBUG_ = kwargs.get('DEBUG', 0)
        self.algorithm = kwargs.get('algorithm', 'follow')
        self.boundary = kwargs.get('markNeighbors', 'search')
        self.output = kwargs.get('output', "")

        if self.output is None:
            self.output = ''
        
        ''' Check to see that all of the meshes exists and that they are
        netcdfs ''' 
        for mesh in mesh_files:
            if os.path.isfile(mesh):
                self.meshes.append(MeshHandler(mesh, 'r', *args, **kwargs))

                if self._DEBUG_ > 0:
                    print("DEBUG: ", mesh, " is a valid NetCDF File\n")
            else:
                print("ERROR: Mesh file was not found", mesh)
                sys.exit(-1)

        ''' Check to see the points file exists and if it exists, then parse it
        and see that is is specified correctly! '''
        if os.path.isfile(region):
            self.region_file = region
        if regionFormat == 'points':
            self.regionSpec = RegionSpec(method=regionFormat, *args, **kwargs)
            self.regionFormat = 'points'
        elif regionFormat == 'shape' .OR. regionFormat == 'shapeFile':
            self.regionSpec = RegionSpec(method=regionFormat, *args, **kwargs)
            self.regionFormat = 'shape'
        else:
            raise NotImplementedError("REGION SPEC IS NOT IMPLMENTED "
                                      "- IMPEMTED IT!")

        ''' Choose the algorithm to choose boundary points '''
        if self.algorithm == 'follow':
            self.mark_boundry = self.follow_the_line

        ''' Choose the algorithm to mark relaxation region '''
        if self.boundary == None:
            # Possibly faster for larger regions
            self.mark_neighbors = self._mark_neighbors
        elif self.boundary == 'search':
            # Possibly faster for smaller regions
            self.mark_neighbors = self._mark_neighbors_search
        

        

    def gen_region(self, *args, **kwargs):
        ''' gen_region

        '''
        # Call the regionSpec to generate `name, in_point, points`
        name, inPoint, points = self.regionSpec.gen_spec(self.region_file)

        if self._DEBUG_ > 0:
            print("DEBUG: Region Spec has been generated")
            print("DEBUG: Name: ", name)
            print("DEBUG: in_point: ", inPoint)
            print("DEBUG: points: ", points)


        # For each mesh, mark the boundary
        for mesh in self.meshes:
            # Mark the boundary cells
            bdyMaskCell, globalBdyCellsIDs, inCell = self.mark_boundry(mesh, 
                                                                       inPoint, 
                                                                       points)
            # Flood fill from the inside point 
            self.flood_fill(mesh, inCell, bdyMaskCell)

            if self._DEBUG_ > 0:
                print("DEBUG: ", globalBdyCellsIDs)

            # Mark the neighbors
            for layer in range(1, self.num_boundary_layers + 1):
                if self._DEBUG_ > 3:
                    print("Debug: Layer: ", layer)

                self.mark_neighbors(mesh, layer, bdyMaskCell, inCell=inCell)

            # Mark the edges
            bdyMaskEdge = self.mark_edges(mesh, bdyMaskCell, *args, **kwargs)

            if self._DEBUG_ > 5:
                print("BdyMaskEdge: ", bdyMaskEdge)
                for edge in bdyMaskEdge:
                    print(edge)

            # Mark the verticies
            bdyMaskVertex = self.mark_vertices(mesh, 
                                               bdyMaskCell, 
                                               *args,
                                               **kwargs)

            if self._DEBUG_ > 5:
                print("BdyMaskVertex: ", bdyMaskVertex)
                for vertex in bdyMaskVertex:
                    print(vertex)


    
            if self._DEBUG_ > 4:
                print("DEBUG: bdyMaskCell: ")
                for cells in bdyMaskCell:
                    print("DEBUG: ", cells)

            # Subset the grid into a new region:
            regionFname = self.create_regional_fname(name, mesh)
            regionalMesh = mesh.subset_fields(regionFname, 
                                              bdyMaskCell,
                                              bdyMaskEdge,
                                              bdyMaskVertex,
                                              *args,
                                              **kwargs)
           
            
            print("Created a regional mesh: ", regionFname)
            mesh.mesh.close()
            regionalMesh.mesh.close()


    ''''''''''''''''''''''''''''''
    ''''''''''''''''''''''''''''''

    def create_regional_fname(self, name, mesh):
        nCells = mesh.mesh.dimensions['nCells'].size
        return os.path.join(self.output, name+'.'+str(nCells)+'.nc')


    # Mark_neighbors_search - Faster for smaller regions ??
    def _mark_neighbors_search(self, mesh, nType, bdyMaskCell, *args, **kwargs):
        inCell = kwargs.get('inCell', None)
        if inCell == None:
            print("ERROR: In cell not found within _mark_neighbors_search")

        nEdgesOnCell = mesh.mesh.variables['nEdgesOnCell'][:]
        cellsOnCell = mesh.mesh.variables['cellsOnCell'][:,:]

        stack = [inCell]
        while len(stack) > 0:
            iCell = stack.pop()
            for i in range(nEdgesOnCell[iCell]):
                j = cellsOnCell[iCell, i] - 1
                if nType > bdyMaskCell[j] >= self.INSIDE:
                    bdyMaskCell[j] = -bdyMaskCell[j]
                    stack.append(j)
                elif bdyMaskCell[j] == 0:
                    bdyMaskCell[j] = nType

        bdyMaskCell[:] = abs(bdyMaskCell[:])


    # mark_neighbors - Faster for larger regions ??
    def _mark_neighbors(self, mesh, nType, bdyMaskCell, *args, **kwargs):
        nCells = len(bdyMaskCell)
        nEdgesOnCell = mesh.mesh.variables['nEdgesOnCell'][:]
        cellsOnCell = mesh.mesh.variables['cellsOnCell'][:,:]

        for iCell in range(nCells):
            if bdyMaskCell[iCell] == self.UNMARKED:
                for i in range(nEdgesOnCell[iCell]):
                    v = cellsOnCell[iCell, i] - 1
                    if bdyMaskCell[v] == 0:
                        bdyMaskCell[v] == nType


    def flood_fill(self, mesh, inCell, bdyMaskCell):
        nEdgesOnCell = mesh.mesh.variables['nEdgesOnCell'][:]
        cellsOnCell = mesh.mesh.variables['cellsOnCell'][:,:]

        stack = [inCell]
        while len(stack) > 0:
            iCell = stack.pop()
            for i in range(nEdgesOnCell[iCell]):
                j = cellsOnCell[iCell, i] - 1
                if bdyMaskCell[j] == self.UNMARKED:
                    bdyMaskCell[j] = 1
                    stack.append(j)

    
    def follow_the_line(self, mesh, inPoint, points, *args, **kwargs):
        ''' Mark the nearest cell to each of the cords in points
        as a boundary cell.
        '''

        if self._DEBUG_ > 0: 
            print("DEBUG: Follow the line: ", self._DEBUG_)

        boundaryCells = []
        nCells = mesh.mesh.dimensions['nCells'].size
        latCell = mesh.mesh.variables['latCell'][:]
        lonCell = mesh.mesh.variables['lonCell'][:]
        cellsOnCell = mesh.mesh.variables['cellsOnCell'][:,:]
        maxEdges = mesh.mesh.dimensions['maxEdges'].size
        nEdgesOnCell = mesh.mesh.variables['nEdgesOnCell'][:]
        indexToCellID = mesh.mesh.variables['indexToCellID'][:]
        sphere_radius = mesh.mesh.sphere_radius

        # Find the nearest cells to the list of given boundary points
        for i in range(0, len(points), 2):
            boundaryCells.append(mesh.nearest_cell(points[i],
                                                   points[i+1]))

        # Find the nearest cell to the inside point
        inCell = mesh.nearest_cell(inPoint[0], inPoint[1])

        if self._DEBUG_ > 0:
            print("DEBUG: Boundary Cells: ", boundaryCells)
            print("DEBUG: Inside point: ", inCell)
            print("DEBUG: Sphere Radius: ", sphere_radius)


        # Create the bdyMask fields
        # TODO: Update this with dtype,
        # TODO: Update this with order, C or F order in memory?
        bdyMaskCell = np.full(nCells, self.UNMARKED)

        # TODO: 
        # bdyMaskEdge = np.full(nCells, self.UNMARKED)
        # bdyMaskVertex = np.full(nCells, self.UNMARKED)
        
        # Mark the boundary cells that were given as input
        for bCells in boundaryCells:
            bdyMaskCell[bCells] = self.INSIDE

        for i in range(len(boundaryCells)):
            sourceCell = boundaryCells[i]
            targetCell = boundaryCells[(i + 1) % len(boundaryCells)]


            if self._DEBUG_ > 0:
                print("DEBUG: sourceCell: ", sourceCell)
                print("DEBUG: targetCell: ", targetCell)
            

            pta = latlon_to_xyz(latCell[sourceCell], 
                                lonCell[sourceCell], 
                                sphere_radius)
            ptb = latlon_to_xyz(latCell[targetCell],
                                lonCell[targetCell],
                                sphere_radius)
        
            pta = np.cross(pta, ptb)
            temp = np.linalg.norm(pta)
            cross = pta / temp
            iCell = sourceCell
            while iCell != targetCell:
                print('iCell: ', iCell)
                bdyMaskCell[iCell] = self.INSIDE
                minangle = np.Infinity
                mindist = sphere_distance(latCell[iCell], 
                                          lonCell[iCell],
                                          latCell[targetCell],
                                          lonCell[targetCell],
                                          sphere_radius)
                for j in range(nEdgesOnCell[iCell]):
                    v = cellsOnCell[iCell,j] - 1
                    dist = sphere_distance(latCell[v],
                                           lonCell[v],
                                           latCell[targetCell],
                                           lonCell[targetCell],
                                           sphere_radius)
                    if dist > mindist:
                        continue
                    pt = latlon_to_xyz(latCell[v], lonCell[v], sphere_radius)
                    angle = np.dot(pta, pt)
                    angle = abs(0.5 * np.pi - np.arccos(angle))
                    if angle < minangle:
                        minangle = angle
                        k = v
                iCell = k

        if self._DEBUG_ > 2:
            print("DEBUG: bdyMaskCells: ", bdyMaskCell)

        return (bdyMaskCell, 
                indexToCellID[np.where(bdyMaskCell != 0)], 
                inCell)

    def mark_edges(self, mesh, bdyMaskCell, *args, **kwargs):

        nCells = mesh.mesh.dimensions['nCells'].size
        nEdges = mesh.mesh.dimensions['nEdges'].size
        cellsOnEdge = mesh.mesh.variables['cellsOnEdge'][:,:]
        bdyMaskEdge = np.full(nEdges, self.UNMARKED)
        
        for edge in range(nEdges):
            if bdyMaskCell[cellsOnEdge[edge, 0]] == self.UNMARKED:
                bdyMaskEdge[edge] = bdyMaskCell[cellsOnEdge[edge, 1] - 1]
            elif bdyMaskEdge[cellsOnEdge[edge, 1]] == self.UNMARKED:
                bdyMaskEdge[edge] = bdyMaskCell[cellsOnEdge[edge, 0]]
            else:
                bdyMaskEdge[edge] = min(bdyMaskCell[cellsOnEdge[edge,0]],
                                        bdyMaskCell[cellsOnEdge[edge,1]])

        return bdyMaskEdge


    def mark_vertices(self, mesh, bdyMaskCell, *args, **kwargs):

        nVertices = mesh.mesh.dimensions['nVertices'].size
        vertexDegree = mesh.mesh.dimensions['vertexDegree'].size
        cellsOnVertex = mesh.mesh.variables['cellsOnVertex'][:,:]

        bdyMaskVertex = np.full(nVertices, self.UNMARKED)

        for vertex in range(nVertices):
            for deg in range(vertexDegree):
                if bdyMaskCell[cellsOnVertex[vertex, deg] - 1] == self.UNMARKED:
                    continue
                elif bdyMaskVertex[vertex] == self.UNMARKED:
                    bdyMaskVertex[vertex] = bdyMaskCell[cellsOnVertex[vertex,
                                                                      deg]]

                if (  bdyMaskCell[cellsOnVertex[vertex, deg]] 
                    < bdyMaskVertex[vertex]):
                    
                    bdyMaskVertex[vertex] = bdyMaskCell[cellsOnVertex[vertex,
                                                                      deg]]
        return bdyMaskVertex
