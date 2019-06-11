from __future__ import absolute_import, division, print_function
import argparse
import os
import sys

import numpy as np

from limited_area.mesh import MeshHandler
from limited_area.mesh import latlon_to_xyz
from limited_area.mesh import sphere_distance
from limited_area.region_spec import RegionSpec

class LimitedArea():
    """ Facilitate creating a regional MPAS mesh from a global MPAS mesh  """
    num_boundary_layers = 8
    INSIDE = 1
    UNMARKED = 0

    def __init__(self,
                 mesh_files,
                 region,
                 regionFormat='points',
                 *args,
                 **kwargs):
        """ Init function for Limited Area

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
        DEBUG         -- Debug value used to turn on debug output, default == 0
        output        -- Optional name to append to regional mesh filename
        markNeighbors -- Algorithm choice for choosing relaxation layers - Default
                         is mark neighbor serach
        """ 
        self.meshes = []

        # Keyword arguments
        self._DEBUG_ = kwargs.get('DEBUG', 0)
        self.boundary = kwargs.get('markNeighbors', 'search')
        self.output = kwargs.get('output', "")

        if self.output is None:
            self.output = ''

        # Check to see that all of the meshes exists and that they are
        # valid netCDF files.
        for mesh in mesh_files:
            if os.path.isfile(mesh):
                self.meshes.append(MeshHandler(mesh, 'r', *args, **kwargs))
            else:
                print("ERROR: Mesh file was not found", mesh)
                sys.exit(-1)

        # Check to see the points file exists and if it exists, then parse it
        # and see that is is specified correctly!
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

        # Choose the algorithm to mark relaxation region
        if self.boundary == None:
            # Possibly faster for larger regions
            self.mark_neighbors = self._mark_neighbors
        elif self.boundary == 'search':
            # Possibly faster for smaller regions
            self.mark_neighbors = self._mark_neighbors_search
        
        
    def gen_region(self, *args, **kwargs):
        """ Generate the boundary region of the given region for the given mesh(es). """

        # Call the regionSpec to generate `name, in_point, points`
        name, inPoint, points = self.regionSpec.gen_spec(self.region_file, **kwargs)

        if self._DEBUG_ > 0:
            print("DEBUG: Region Spec has been generated")
            print("DEBUG: Region Name: ", name)
            print("DEBUG: In Point: ", inPoint)
            print("DEBUG: # of points: ", len(points))

        # For each mesh, create a regional mesh and save it
        for mesh in self.meshes:
            print('\n')
            print('Creating a regional mesh of ', mesh.fname)

            # Mark the boundary cells
            print('Marking boundary cells ...')
            bdyMaskCell, globalBdyCellsIDs, inCell = self.mark_boundary(mesh,
                                                                        inPoint,
                                                                        points)
            # Flood fill from the inside point 
            print('Filling region ...')
            bdyMaskCell = self.flood_fill(mesh, inCell, bdyMaskCell)

            # Mark the neighbors
            print('Creating boundary laryer:', end=' '); sys.stdout.flush()
            for layer in range(1, self.num_boundary_layers + 1):
                print(layer, ' ...', end=' '); sys.stdout.flush()
                self.mark_neighbors(mesh, layer, bdyMaskCell, inCell=inCell)
            print('DONE!')

            if self._DEBUG_ > 2:
                print("DEBUG: bdyMaskCells count:")
                print("DEBUG: 0: ", len(bdyMaskCell[bdyMaskCell == 0]))
                print("DEBUG: 1: ", len(bdyMaskCell[bdyMaskCell == 1]))
                print("DEBUG: 2: ", len(bdyMaskCell[bdyMaskCell == 2]))
                print("DEBUG: 3: ", len(bdyMaskCell[bdyMaskCell == 3]))
                print("DEBUG: 4: ", len(bdyMaskCell[bdyMaskCell == 4]))
                print("DEBUG: 5: ", len(bdyMaskCell[bdyMaskCell == 5]))
                print("DEBUG: 6: ", len(bdyMaskCell[bdyMaskCell == 6]))
                print("DEBUG: 7: ", len(bdyMaskCell[bdyMaskCell == 7]))
                print("DEBUG: 8: ", len(bdyMaskCell[bdyMaskCell == 8]))
                print('\n')

            bdyMaskCell_cp = bdyMaskCell

            # Mark the edges
            print('Markin region edges ...')
            bdyMaskEdge = self.mark_edges(mesh, 
                                          bdyMaskCell, 
                                          *args, 
                                          **kwargs)

            # Mark the verticies
            print('Markin region verteices...')
            bdyMaskVertex = self.mark_vertices(mesh, 
                                               bdyMaskCell, 
                                               *args,
                                               **kwargs)


            # Subset the grid into a new region:
            print('Subseting mesh fields into the specified region mesh...')
            regionFname = self.create_regional_fname(name, mesh, output=self.output)
            regionalMesh = mesh.subset_fields(regionFname,
                                              bdyMaskCell,
                                              bdyMaskEdge,
                                              bdyMaskVertex,
                                              inside=self.INSIDE,
                                              unmarked=self.UNMARKED,
                                              *args,
                                              **kwargs)

            print('Copying global attributes...')
            mesh.copy_global_attributes(regionalMesh)

            print("Created a regional mesh: ", regionFname)

            print('Creating graph partition file...', end=' '); sys.stdout.flush()
            regionalMesh.create_graph_file(self.create_partiton_fname(name, 
                                                                      mesh, 
                                                                      output=self.output))

            mesh.mesh.close()
            regionalMesh.mesh.close()

    def create_partiton_fname(self, name, mesh, **kwargs):
        """ Generate the filename for the regional graph.info file"""
        output = kwargs.get('output', None)

        if output:
            return output+'.graph.info'
        else:
            nCells = mesh.mesh.dimensions['nCells'].size
            return name+'.'+str(nCells)+'.graph.info'
        

    def create_regional_fname(self, name, mesh, **kwargs):
        """ Generate the filename for the regional mesh file """
        output = kwargs.get('output', None)

        if output:
            return output 
        else:
            nCells = mesh.mesh.dimensions['nCells'].size
            return name+'.'+str(nCells)+'.grid.nc'


    # Mark_neighbors_search - Faster for smaller regions ??
    def _mark_neighbors_search(self, mesh, layer, bdyMaskCell, *args, **kwargs):
        """ Mark the relaxation layers using a search and return an updated bdyMaskCell with
        those relaxation layers
        
        mesh        -- The global MPAS mesh
        layer       -- The relaxation layer
        bdyMaskCell -- The global mask marking the regional cell subset
        inCell      -- A point that is inside the regional area

        """
        inCell = kwargs.get('inCell', None)
        if inCell == None:
            print("ERROR: In cell not found within _mark_neighbors_search")

        nEdgesOnCell = mesh.mesh.variables['nEdgesOnCell'][:]
        cellsOnCell = mesh.mesh.variables['cellsOnCell'][:, :]

        stack = [inCell]
        while len(stack) > 0:
            iCell = stack.pop()
            for i in range(nEdgesOnCell[iCell]):
                j = cellsOnCell[iCell, i] - 1
                if layer > bdyMaskCell[j] >= self.INSIDE:
                    bdyMaskCell[j] = -bdyMaskCell[j]
                    stack.append(j)
                elif bdyMaskCell[j] == 0:
                    bdyMaskCell[j] = layer 

        bdyMaskCell[:] = abs(bdyMaskCell[:])


    # mark_neighbors - Faster for larger regions ??
    def _mark_neighbors(self, mesh, nType, bdyMaskCell, *args, **kwargs):
        """ """
        nCells = len(bdyMaskCell)
        nEdgesOnCell = mesh.mesh.variables['nEdgesOnCell'][:]
        cellsOnCell = mesh.mesh.variables['cellsOnCell'][:, :]

        for iCell in range(nCells):
            if bdyMaskCell[iCell] == self.UNMARKED:
                for i in range(nEdgesOnCell[iCell]):
                    v = cellsOnCell[iCell, i] - 1
                    if bdyMaskCell[v] == 0:
                        bdyMaskCell[v] == nType


    def flood_fill(self, mesh, inCell, bdyMaskCell):
        """ Mark the interior points of the regional mesh and return and updated
        bdyMaskCell.

        mesh        -- Global MPAS Mesh
        inCell      -- A point that is inside the specified region
        bdyMaskCell -- The global mask marking which global cells are interior, relaxation
                       and those that are outside.
        """
        if self._DEBUG_ > 1:
            print("DEBUG: Flood filling with flood_fill!")
        nEdgesOnCell = mesh.mesh.variables['nEdgesOnCell'][:]
        cellsOnCell = mesh.mesh.variables['cellsOnCell'][:, :]

        stack = [inCell]
        while len(stack) > 0:
            iCell = stack.pop()
            for i in range(nEdgesOnCell[iCell]):
                j = cellsOnCell[iCell, i] - 1
                if bdyMaskCell[j] == self.UNMARKED:
                    bdyMaskCell[j] = self.INSIDE
                    stack.append(j)

        return bdyMaskCell


    def mark_edges(self, mesh, bdyMaskCell, *args, **kwargs):
        """ Mark the edges that are in the specified region and return
        bdyMaskEdge.
        
        mesh -
        bdyMaskCell -
        """
        nEdges = mesh.mesh.dimensions['nEdges'].size
        cellsOnEdge = mesh.mesh.variables['cellsOnEdge'][:]
        bdyMaskEdge = np.zeros(nEdges, dtype=np.dtype('i'))
        np.set_printoptions(threshold=np.inf)

        for i in range(nEdges):
            cells = cellsOnEdge[i,:]
            if np.all(bdyMaskCell[cells - 1 ] == 0):
                bdyMaskEdge[i] = 0
            elif np.any(bdyMaskCell[cells - 1] == 0):
                cellMasks = bdyMaskCell[cells - 1]
                bdyMaskEdge[i] = np.min(cellMasks[cellMasks > 0])
            else:
                bdyMaskEdge[i] = np.min(bdyMaskCell[cells - 1])

        if self._DEBUG_ > 2:
            print("DEBUG: bdyMaskEdges count:")
            print("DEBUG: 0: ", len(bdyMaskEdge[bdyMaskEdge == 0]))
            print("DEBUG: 1: ", len(bdyMaskEdge[bdyMaskEdge == 1]))
            print("DEBUG: 2: ", len(bdyMaskEdge[bdyMaskEdge == 2]))
            print("DEBUG: 3: ", len(bdyMaskEdge[bdyMaskEdge == 3]))
            print("DEBUG: 4: ", len(bdyMaskEdge[bdyMaskEdge == 4]))
            print("DEBUG: 5: ", len(bdyMaskEdge[bdyMaskEdge == 5]))
            print("DEBUG: 6: ", len(bdyMaskEdge[bdyMaskEdge == 6]))
            print("DEBUG: 7: ", len(bdyMaskEdge[bdyMaskEdge == 7]))
            print("DEBUG: 8: ", len(bdyMaskEdge[bdyMaskEdge == 8]))
            print('\n')

        return bdyMaskEdge


    def mark_vertices(self, mesh, bdyMaskCell, *args, **kwargs):
        """ Mark the vertices that are in the spefied region and return
        bdyMaskVertex."""
        nVertices = mesh.mesh.dimensions['nVertices'].size
        vDegree = mesh.mesh.dimensions['vertexDegree'].size
        cellsOnVertex = mesh.mesh.variables['cellsOnVertex'][:]
        bdyMaskVertex = np.zeros(nVertices, dtype=np.dtype('i'))

        for i in range(nVertices):
            cells = cellsOnVertex[i,:]
            if np.all(bdyMaskCell[cells - 1] == 0):
                bdyMaskVertex[i] = 0
            elif np.any(bdyMaskCell[cells - 1] == 0):
                cellMasks = bdyMaskCell[cells - 1]
                bdyMaskVertex[i] = np.min(cellMasks[cellMasks > 0])
            else:
                bdyMaskVertex[i] = np.min(bdyMaskCell[cells - 1])

        if self._DEBUG_ > 2:
            print("DEBUG: bdyMaskVertex count:")
            print("DEBUG: 0: ", len(bdyMaskVertex[bdyMaskVertex == 0]))
            print("DEBUG: 1: ", len(bdyMaskVertex[bdyMaskVertex == 1]))
            print("DEBUG: 2: ", len(bdyMaskVertex[bdyMaskVertex == 2]))
            print("DEBUG: 3: ", len(bdyMaskVertex[bdyMaskVertex == 3]))
            print("DEBUG: 4: ", len(bdyMaskVertex[bdyMaskVertex == 4]))
            print("DEBUG: 5: ", len(bdyMaskVertex[bdyMaskVertex == 5]))
            print("DEBUG: 6: ", len(bdyMaskVertex[bdyMaskVertex == 6]))
            print("DEBUG: 7: ", len(bdyMaskVertex[bdyMaskVertex == 7]))
            print("DEBUG: 8: ", len(bdyMaskVertex[bdyMaskVertex == 8]))
            print('\n')

        return bdyMaskVertex
    

    # Mark Boundary points
    def mark_boundary(self, mesh, inPoint, points, *args, **kwargs):
        """ Mark the nearest cell to each of the cords in points
        as a boundary cell.

        mesh -
        inPoint -
        points -

        """
        if self._DEBUG_ > 0: 
            print("DEBUG: Marking the boundary points: ")

        boundaryCells = []
        nCells = mesh.mesh.dimensions['nCells'].size
        latCell = mesh.mesh.variables['latCell'][:]
        lonCell = mesh.mesh.variables['lonCell'][:]
        cellsOnCell = mesh.mesh.variables['cellsOnCell'][:, :]
        maxEdges = mesh.mesh.dimensions['maxEdges'].size
        nEdgesOnCell = mesh.mesh.variables['nEdgesOnCell'][:]
        indexToCellID = mesh.mesh.variables['indexToCellID'][:]
        sphere_radius = mesh.mesh.sphere_radius

        # Find the nearest cells to the list of given boundary points
        for i in range(0, len(points), 2):
            boundaryCells.append(mesh.nearest_cell(points[i],
                                                   points[i + 1]))


        # Find the nearest cell to the inside point
        inCell = mesh.nearest_cell(inPoint[0], inPoint[1])

        if self._DEBUG_ > 0:
            print("DEBUG: Num Boundary Cells: ", len(boundaryCells))
            print("DEBUG: Inside Cell: ", inCell)

        # Create the bdyMask fields
        bdyMaskCell = np.full(nCells, self.UNMARKED)

        # Mark the boundary cells that were given as input
        for bCells in boundaryCells:
            bdyMaskCell[bCells] = self.INSIDE

        # For each boundaryCells, mark the current cell as the source cell
        # and the next (or the first element if the current is the last) as 
        # the target cell.
        #
        # Then, determine the great-arc angle between the source and taget
        # cell, and then for each cell, starting at the source cell, 
        # calculate the great-arc angle between the cells on the current
        # cell and the target cell, and then add the cell with the smallest
        # angle.
        for i in range(len(boundaryCells)):
            sourceCell = boundaryCells[i]
            targetCell = boundaryCells[(i + 1) % len(boundaryCells)]

            # If we are already at the next target cell, there is no need
            # to connect sourceCell with targetCell, and we can skip to
            # the next pair of boundary points
            if sourceCell == targetCell:
                continue

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
                bdyMaskCell[iCell] = self.INSIDE
                minangle = np.Infinity
                mindist = sphere_distance(latCell[iCell], 
                                          lonCell[iCell],
                                          latCell[targetCell],
                                          lonCell[targetCell],
                                          sphere_radius)
                for j in range(nEdgesOnCell[iCell]):
                    v = cellsOnCell[iCell, j] - 1
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

        return (bdyMaskCell, 
                indexToCellID[np.where(bdyMaskCell != self.UNMARKED)], 
                inCell)

