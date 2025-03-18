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
                 mesh_file,
                 region,
                 plotting,
                 regionFormat='points',
                 format='NETCDF3_64BIT_OFFSET',
                 *args,
                 **kwargs):
        """ Init function for Limited Area

        Check to see if mesh file exists and it is the correct type. Check to
        see that the region file exist and finally set the regionSpec to the
        requested regionFormat


        Keyword arguments:
        mesh_files   -- Path to a valid MPAS Mesh file
        region       -- Path to pts file region specification 

        DEBUG         -- Debug value used to turn on debug output, default == 0
        markNeighbors -- Algorithm choice for choosing relaxation layers - Default
                         is mark neighbor search
        """ 

        # Keyword arguments
        self._DEBUG_ = kwargs.get('DEBUG', 0)
        self.boundary = kwargs.get('markNeighbors', 'search')
        self.plotting = plotting
        self.cdf_format = format

        # Check to see that all of the meshes exists and that they are
        # valid netCDF files.
        if not self.plotting:
            if os.path.isfile(mesh_file):
                self.mesh = MeshHandler(mesh_file, 'r', *args, **kwargs)
            else:
                print("ERROR: Mesh file was not found", mesh_file)
                sys.exit(-1)

        # Check to see the points file exists and if it exists, then parse it
        # and see that is is specified correctly!
        self.region_file = region
        self.regionSpec = RegionSpec(*args, **kwargs)

        # Choose the algorithm to mark relaxation region
        if self.boundary == None:
            # Possibly faster for larger regions
            self.mark_neighbors = self._mark_neighbors
        elif self.boundary == 'search':
            # Possibly faster for smaller regions
            self.mark_neighbors = self._mark_neighbors_search
        
        
    def gen_region(self, *args, **kwargs):
        """ Generate the boundary region of the given region for the given mesh(es). """

        # Call the regionSpec to generate `name, in_point, boundaries`
        name, inPoint, boundaries= self.regionSpec.gen_spec(self.region_file, self.plotting, **kwargs)

        if self.plotting:
            self.plot_region(inPoint, boundaries, name)
            return

        if self._DEBUG_ > 0:
            print("DEBUG: Region Spec has been generated")
            print("DEBUG: Region Name: ", name)
            print("DEBUG: In Point: ", inPoint)
            print("DEBUG: # of boundaries: ", len(boundaries))

        # For each mesh, create a regional mesh and save it
        print('\n')
        print('Creating a regional mesh of ', self.mesh.fname)

        # Mark boundaries
        # A specification may have multiple, discontiguous boundaries,
        # so, create a unmarked, filled bdyMaskCell and pass it to
        # mark_boundary for each boundary.
        print('Marking ', end=''); sys.stdout.flush()
        bdyMaskCell = np.full(self.mesh.nCells, self.UNMARKED)
        i = 1
        for boundary in boundaries:
            print("boundary ", i, "... ", end=''); sys.stdout.flush(); i += 1
            bdyMaskCell = self.mark_boundary(self.mesh, boundary, bdyMaskCell)

        # Find the nearest cell to the inside point
        inCell = self.mesh.nearest_cell(inPoint[0], inPoint[1])

        # Flood fill from the inside point
        print('\nFilling region ...')
        bdyMaskCell = self.flood_fill(self.mesh, inCell, bdyMaskCell)

        # Mark the neighbors
        print('Creating boundary layer:', end=' '); sys.stdout.flush()
        for layer in range(1, self.num_boundary_layers + 1):
            print(layer, ' ...', end=' '); sys.stdout.flush()
            self.mark_neighbors(self.mesh, layer, bdyMaskCell, inCell=inCell)
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
        print('Marking region edges ...')
        bdyMaskEdge = self.mark_edges(self.mesh,
                                      bdyMaskCell,
                                      *args,
                                      **kwargs)

        # Mark the vertices
        print('Marking region vertices...')
        bdyMaskVertex = self.mark_vertices(self.mesh,
                                           bdyMaskCell,
                                           *args,
                                           **kwargs)


        # Subset the grid into a new region:
        print('Subsetting mesh fields into the specified region mesh...')
        regionFname = self.create_regional_fname(name, self.mesh)
        regionalMesh = self.mesh.subset_fields(regionFname,
                                          bdyMaskCell,
                                          bdyMaskEdge,
                                          bdyMaskVertex,
                                          inside=self.INSIDE,
                                          unmarked=self.UNMARKED,
                                          format=self.cdf_format,
                                          *args,
                                          **kwargs)

        print('Copying global attributes...')
        self.mesh.copy_global_attributes(regionalMesh)

        print("Created a regional mesh: ", regionFname)

        print('Creating graph partition file...', end=' '); sys.stdout.flush()
        graphFname = regionalMesh.create_graph_file(self.create_partiton_fname(name, self.mesh,))
        print(graphFname)

        self.mesh.mesh.close()
        regionalMesh.mesh.close()

        return regionFname, graphFname

    def create_partiton_fname(self, name, mesh, **kwargs):
        """ Generate the filename for the regional graph.info file"""
        return name+'.graph.info'
        

    def create_regional_fname(self, name, mesh, **kwargs):
        """ Generate the filename for the regional mesh file. Depending on what "type" of MPAS file
        we are using, either static, grid or init, try to rename the region file as that type (i.e.
        x1.2562.static.nc becomes name.static.nc).
        
        If a file name is ambiguous, or the file name does not contain: static, init, or grid,
        rename the region file to be region. """
        # Static files
        if 'static' in mesh.fname and not ('grid' in mesh.fname or 'init' in mesh.fname):
            meshType = 'static'
        # Grid files
        elif 'grid' in mesh.fname and not ('static' in mesh.fname or 'init' in mesh.fname):
            meshType = 'grid'
        # Initialization Data
        elif 'init' in mesh.fname and not ('static' in mesh.fname or 'grid' in mesh.fname):
            meshType = 'init'
        else:
            meshType = 'region'

        return name+'.'+meshType+'.nc'


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

        stack = [inCell]
        while len(stack) > 0:
            iCell = stack.pop()
            for i in range(mesh.nEdgesOnCell[iCell]):
                j = mesh.cellsOnCell[iCell, i] - 1
                if layer > bdyMaskCell[j] >= self.INSIDE:
                    bdyMaskCell[j] = -bdyMaskCell[j]
                    stack.append(j)
                elif bdyMaskCell[j] == 0:
                    bdyMaskCell[j] = layer 

        bdyMaskCell[:] = abs(bdyMaskCell[:])


    # mark_neighbors - Faster for larger regions ??
    def _mark_neighbors(self, mesh, nType, bdyMaskCell, *args, **kwargs):
        """ Mark a relaxation layers of nType

        mesh        -- The global MPAS mesh
        nType       -- The current relaxation cell that will be marked on bdyMaskCell
        bdyMaskCell -- The global mask marking the regional cell subset
        """

        for iCell in range(mesh.nCells):
            if bdyMaskCell[iCell] == self.UNMARKED:
                for i in range(mesh.nEdgesOnCell[iCell]):
                    v = mesh.cellsOnCell[iCell, i] - 1
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
            print("DEBUG: Flood filling with flood_fill")

        stack = [inCell]
        while len(stack) > 0:
            iCell = stack.pop()
            for i in range(mesh.nEdgesOnCell[iCell]):
                j = mesh.cellsOnCell[iCell, i] - 1
                if bdyMaskCell[j] == self.UNMARKED:
                    bdyMaskCell[j] = self.INSIDE
                    stack.append(j)

        return bdyMaskCell


    def mark_edges(self, mesh, bdyMaskCell, *args, **kwargs):
        """ Mark the edges that are in the specified region and return
        bdyMaskEdge. """

        bdyMaskEdge = bdyMaskCell[mesh.cellsOnEdge[:,:]-1].min(axis=1)
        bdyMaskEdge = np.where(bdyMaskEdge > 0,
                               bdyMaskEdge,
                               bdyMaskCell[mesh.cellsOnEdge[:,:]-1].max(axis=1))

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

        bdyMaskVertex = bdyMaskCell[mesh.cellsOnVertex[:,:]-1].min(axis=1)
        bdyMaskVertex = np.where(bdyMaskVertex > 0,
                                 bdyMaskVertex,
                                 bdyMaskCell[mesh.cellsOnVertex[:,:]-1].max(axis=1))

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
    def mark_boundary(self, mesh, points, bdyMaskCell, *args, **kwargs):
        """ Mark the nearest cell to each of the cords in points
        as a boundary cell and return bdyMaskCell.

        mesh - The global mesh
        inPoint - A point that lies within the regional area
        points - A list of points that define the boundary of the desired
                 region as flatten list of lat, lon coordinates. i.e:

                 [lat0, lon0, lat1, lon1, lat2, lon2, ..., latN, lonN]
        """
        if self._DEBUG_ > 0:
            print("DEBUG: Marking the boundary points: ")

        boundaryCells = []

        # Find the nearest cells to the list of given boundary points
        for i in range(0, len(points), 2):
            boundaryCells.append(mesh.nearest_cell(points[i],
                                                   points[i + 1]))



        if self._DEBUG_ > 0:
            print("DEBUG: Num Boundary Cells: ", len(boundaryCells))

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

            pta = latlon_to_xyz(mesh.latCells[sourceCell],
                                mesh.lonCells[sourceCell],
                                1.0)
            ptb = latlon_to_xyz(mesh.latCells[targetCell],
                                mesh.lonCells[targetCell],
                                1.0)
        
            pta = np.cross(pta, ptb)
            temp = np.linalg.norm(pta)
            pta = pta / temp
            iCell = sourceCell
            while iCell != targetCell:
                bdyMaskCell[iCell] = self.INSIDE
                minangle = np.Infinity
                mindist = sphere_distance(mesh.latCells[iCell],
                                          mesh.lonCells[iCell],
                                          mesh.latCells[targetCell],
                                          mesh.lonCells[targetCell],
                                          1.0)
                for j in range(mesh.nEdgesOnCell[iCell]):
                    v = mesh.cellsOnCell[iCell, j] - 1
                    dist = sphere_distance(mesh.latCells[v],
                                           mesh.lonCells[v],
                                           mesh.latCells[targetCell],
                                           mesh.lonCells[targetCell],
                                           1.0)
                    if dist > mindist:
                        continue
                    pt = latlon_to_xyz(mesh.latCells[v], mesh.lonCells[v], 1.0)
                    angle = np.dot(pta, pt)
                    angle = abs(0.5 * np.pi - np.arccos(angle))
                    if angle < minangle:
                        minangle = angle
                        k = v
                iCell = k

        return bdyMaskCell


    # Plot the specified region
    def plot_region(self, inPoint, boundaries, region_name):
        """ Create an plot showing the user-specified region.
        The region boundary is plotted as a blue line on an orthographic
        projection with a simple color-filled map background.

        inPoint - A point that lies within the regional area.
        boundaries - A list of lists of latitude and longitude coordinates
                 of a boundary.

                 [[lat0, lon0, lat1, lon1, lat2, lon2, ..., latN, lonN], [...]]

        Note: Coordinates are in radians.
        """

        try:
            import cartopy.crs as ccrs
            import cartopy.feature as cfeature
        except:
            print('Could not import cartopy!')
            return

        try:
            import matplotlib.pyplot as plt
            import matplotlib.ticker as mticker
            import matplotlib.patches as mpatches
        except:
            print('Could not import matplotlib!')
            return

        import math

        rad2deg = 180.0 / math.pi

        if len(boundaries) == 2:
            proj = ccrs.PlateCarree(rad2deg * inPoint[1])
            ax = plt.axes(projection=proj)

            latbdy = np.concatenate((np.flip(boundaries[0][0::2]), boundaries[1][0::2])) * rad2deg
            lonbdy = np.concatenate((np.flip(boundaries[0][1::2]), boundaries[1][1::2])) * rad2deg

            ax.set_extent((-180.0, 180.0, -90.0, 90.0), crs=ccrs.PlateCarree())

        else:
            proj = ccrs.Stereographic(rad2deg * inPoint[0], rad2deg * inPoint[1])
            ax = plt.axes(projection=proj)

            lonbdy, latbdy = self.ccw_order(inPoint[1] * rad2deg, inPoint[0] * rad2deg,
                boundaries[0][1::2] * rad2deg,boundaries[0][0::2] * rad2deg,
                proj)

            # Logic to determine the "extent" (roughly, the diameter in meters)
            # of the region given the latbdy and lonbdy arrays
            # Convert all points from lat,lon coordinates to x,y
            # distances (in meters) from the inPoint
            xy=proj.transform_points(ccrs.PlateCarree(), lonbdy, latbdy)

            # Find the extents of the box bounding these points and the
            # diameter (length of diagonal)
            min_x = min(xy[:,0])
            max_x = max(xy[:,0])
            min_y = min(xy[:,1])
            max_y = max(xy[:,1])
            diameter = math.sqrt((max_y-min_y)**2+(max_x-min_x)**2)

            extent = diameter
            scaling = 0.5 * 1.25
            ax.set_extent([-scaling * extent, scaling * extent, -scaling * extent, scaling * extent], crs=proj)

        # Place latbdy and lonbdy in clockwise order to cause the filled polygon
        # to be everything except the region
        latbdy = latbdy[::-1]
        lonbdy = lonbdy[::-1]

        ax.add_feature(cfeature.OCEAN)
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.1)
        ax.add_feature(cfeature.BORDERS, linewidth=0.1)
        ax.add_feature(cfeature.LAKES, linewidth=0.1)
        ax.add_feature(cfeature.RIVERS, linewidth=0.1)
        ax.add_feature(cfeature.STATES, linewidth=0.1)

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                          linewidth=0.1, color='black', alpha=1.0, linestyle='--')

        xticks = np.arange(-180, 180, 10)
        yticks = np.arange(-90, 90, 10)

        gl.ylocator = mticker.FixedLocator(yticks)
        gl.xlocator = mticker.FixedLocator(xticks)

        # Hack for plotting polygons when all latitudes are the same.
        if np.min(latbdy) == np.max(latbdy):
            if np.min(latbdy) == 0.0:
               latbdy[0] = latbdy[0] + 0.000001
            else:
               latbdy[0] = latbdy[0] * 0.999999

        nbdy = len(lonbdy)
        poly_corners = np.zeros((nbdy+1, 2), np.float64)
        poly_corners[0:nbdy,0] = np.asarray(lonbdy)
        poly_corners[0:nbdy,1] = np.asarray(latbdy)
        poly_corners[nbdy,:] = poly_corners[0,:]

        poly = mpatches.Polygon(poly_corners, closed=True, ec='black', fill=True, lw=0.1, fc='black', alpha=0.3, transform=ccrs.Geodetic())
        ax.add_patch(poly)

        for boundary in boundaries:
            latbdy = boundary[0::2] * rad2deg
            lonbdy = boundary[1::2] * rad2deg

            nbdy = len(lonbdy)

            lon_corners = np.zeros((nbdy+1), np.float64)
            lon_corners[0:nbdy] = np.asarray(lonbdy)
            lon_corners[nbdy] = lon_corners[0]

            lat_corners = np.zeros((nbdy+1), np.float64)
            lat_corners[0:nbdy] = np.asarray(latbdy)
            lat_corners[nbdy] = lat_corners[0]

            plt.plot(lon_corners, lat_corners,
                     color='blue', linewidth=0.5,
                     transform=ccrs.Geodetic()
                     )

        print('Saving plot to region.png')
        plt.savefig('region.png', dpi=150, bbox_inches='tight')


    def plane_angle(self, ax, ay, bx, by, cx, cy):
        """ Computes the angle between the vectors AB and AC.
        """

        import math

        b = [bx - ax, by - ay]
        c = [cx - ax, cy - ay]

        arg = np.dot(b, c) / np.linalg.norm(b) / np.linalg.norm(c)
        arg = max(-1.0,min(1.0, arg))
        theta = math.acos(arg)
        if np.cross([bx - ax, by - ay], [cx - ax, cy - ay]) < 0.0:
            theta = -theta

        return theta


    def ccw_order(self, cenLon, cenLat, lons, lats, proj):
        """ Place arrays of longitude and latitude in counter-clockwise (CCW)
        order. The cenLon and cenLat values define a point that is within the
        polygon described by the lons and lats arrays.

        Latitude and longitude coordinates are in degrees.
        """

        import cartopy.crs as ccrs

        cenLonProj, cenLatProj = proj.transform_point(cenLon, cenLat, ccrs.PlateCarree())
        coordsProj = proj.transform_points(ccrs.PlateCarree(), lons, lats)

        stride = max(1, int(lons.size / 360))

        lonsProj = coordsProj[::stride,0]
        latsProj = coordsProj[::stride,1]

        sumangle = 0.0
        for i in range(lonsProj.size):
            sumangle = sumangle + self.plane_angle(cenLonProj, cenLatProj, lonsProj[i-1], latsProj[i-1], lonsProj[i], latsProj[i])

        if sumangle > 0.0:
            return lons, lats
        else:
            return lons[::-1], lats[::-1]
