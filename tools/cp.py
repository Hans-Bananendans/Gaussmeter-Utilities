"""
Cube Projector - Merged

Merged library version of the Cube Projector framework from
https://github.com/Hans-Bananendans/CubeSat-Solar-Estimator
The reason for merging is for ease of use in other projects.

Code by Johan Monster

This merged file contains the contents of:
cp_vertex.py
cp_face.py
cp_geometry.py
cp_vector.py
cp_frame.py
cp_plotting.py
cp_utilities.py
"""

import numpy as np
from numpy import sin, cos
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mp3d

class Vertex:
    """Custom implentation of a 3D point expressed in cartesians.

    Usage: Vertex has 3D carthesian coordinates, and some knowledge of its
    parent. The parent can be used to guide certain behaviour.

    Construction examples:
        p1 = Vertex()
            Creates a Vertex instance 'p1' at XYZ(0, 0, 0).

        p1 = Vertex([2, 4, 1])
            Creates a Vertex instance 'p1' at XYZ (2, 4, 1).

        p1 = Vertex([1, 1, 1], parenttype='frame')
            Creates a Vertex instance 'p1' at XYZ (1, 1, 1) and manually
            set the vertex parenttype to 'frame' (generally not necessary).

    """

    def __init__(self, xyz=(0, 0, 0), parenttype="global"):

        if not (isinstance(xyz, list) or isinstance(xyz, np.ndarray)):
            raise TypeError("Vertex constructor was not supplied with correct\
                            Vertex coordinates! Use a list or numpy.ndarray!")

        self.x = xyz[0]
        self.y = xyz[1]
        self.z = xyz[2]
        self.parenttype = "global"

        # Re-assign self.parenttype to check for validity of given parenttype
        if parenttype != "global":
            self.set_parenttype(parenttype)

    def set_parenttype(self, new_parenttype):
        """Sets own parenttype to a specified parenttype. First verifies
        specified parenttype.
        """
        if new_parenttype in ["global", "frame", "geometry", "face", "vector"]:
            self.parenttype = new_parenttype
        else:
            raise ValueError("The parenttype cannot be anything other than: \
                             'global','frame','geometry','face','vector'!")

    def translate(self, dx=0., dy=0., dz=0.):
        """Move the vertex' carthesian components by a relative amount.
        """
        self.x += dx
        self.y += dy
        self.z += dz

    def rotate(self, a, b, c, cor=None, seq='321'):
        """ Rotates the vertex around a point in space.
            a, b, c is Euler angle of rotation around local x, y, z axes
                expressed in radians

            seq is a string with the rotation sequence, e.g. '321' for:
                Rz(c).Ry(b).Rx(a).vertex
                Default rotation sequence: '321' (zyx)

            cor is the centre of rotation, which can be specified as
                a Vertex, or a list or numpy.ndarray of local coordinates.
        """

        # ==== COR argument handling ====
        # Case: COR is given as a list/ndarray of coordinates
        if isinstance(cor, list) or isinstance(cor, np.ndarray):
            if len(cor) == 3:
                cor = Vertex(list(cor))

        # Case: COR is given as a Vertex
        elif isinstance(cor, Vertex):
            pass

        # All other cases
        else:
            raise TypeError("Error when rotating vertex. Centre of rotation \
                            argument not provided correctly!")

        # ==== Define relevant rotation matrices ====
        Rx = np.array([[1, 0, 0],
                       [0, cos(a), -sin(a)],
                       [0, sin(a), cos(a)]])

        Ry = np.array([[cos(b), 0, sin(b)],
                       [0, 1, 0],
                       [-sin(b), 0, cos(b)]])

        Rz = np.array([[cos(c), -sin(c), 0],
                       [sin(c), cos(c), 0],
                       [0, 0, 1]])

        # ==== Perform rotation ====
        for axis in reversed(seq):

            # If no Centre of Rotation (cor) is given, take frame origin.
            # If not, subtract coordinates of cor first before applying
            # the rotation, and then reverse this afterwards.
            if cor:
                self.translate(-1 * cor.x, -1 * cor.y, -1 * cor.z)

            if axis == '1':
                (self.x, self.y, self.z) = \
                    np.dot(Rx, np.array([self.x, self.y, self.z]))

            if axis == '2':
                (self.x, self.y, self.z) = \
                    np.dot(Ry, np.array([self.x, self.y, self.z]))

            if axis == '3':
                (self.x, self.y, self.z) = \
                    np.dot(Rz, np.array([self.x, self.y, self.z]))

            # Reverse temporary translation.
            if cor:
                self.translate(1 * cor.x, 1 * cor.y, 1 * cor.z)

    def xyz(self):
        """Returns a numpy array with the local x, y, z coordinates.

            Output: np.array([x, y, z])
        """
        return np.array([self.x, self.y, self.z])

    def readout(self, dec=4):
        """More detailed information about Vertex.
           Use dec to define the number of decimals before rounding.
           Set dec to -1 to disable rounding.
           """
        if dec == -1:  # Set dec to -1 to disable rounding
            localcoords = self.xyz()
        else:
            localcoords = np.round(self.xyz(), dec)
        print("Local coordinates:   {}".format(localcoords))


class Face:
    """Class for a quadrilateral face (4-point polygon).

       A face consists of four vertices that make up the corners of a plane.
       All four vertices must be coplanar, and the class constructor will
       verify this. The order in which the vertices are supplied DOES matter,
       as the faces have one side only (as opposed to a front and a backside)

       For best results the corners must be supplied sequentially and
       counterclockwise, i.e.:

          4  o---o 3
             |   |
          1  o---o 2

       This will create a face that has a front area, and the normal of the
       face will point upward (towards the reader). If supplied in clockwise
       direction instead, this face would have identical points, but its
       normal would point down (away from the reader, into the screen).

       Construction example:

           f1 = Face(Vertex([ 1, 0, 0]),
                     Vertex([ 1, 1, 0]),
                     Vertex([-1, 1, 0]),
                     Vertex([-1,-1, 0])
                     )
               Creates an instance of the Face class named 'f1', which forms
               a square face in the XY-plane, with sides of length 2, centred
               on the local origin.
       """

    def __init__(self, p1: Vertex, p2: Vertex, p3: Vertex, p4: Vertex,
                 parenttype="global"):

        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4

        # Ensure provided vertices are coplanar
        if self.check_coplanarity() is False:
            raise ValueError("Error! Set of four points is not coplanar!")

        self.parenttype = "global"

        # Re-assign self.parenttype to check for validity of given parenttype
        if parenttype != "global":
            self.set_parenttype(parenttype)

        # Let child vertices know their parenttype is 'face'
        for vertex in self.vertices():
            vertex.set_parenttype("face")

    def set_parenttype(self, new_parenttype):
        """Sets own parenttype to a specified parenttype. First verifies
        specified parenttype.
        """
        if new_parenttype in ["global", "frame", "geometry"]:
            self.parenttype = new_parenttype

        else:
            raise ValueError("The parenttype cannot be anything other than: \
                             'global', 'frame', 'geometry'!")

    def check_coplanarity(self):
        """Check if all points are coplanar, by investigating determinant
               of [p12, p13, p14]. Points are coplanar only if it is zero."""
        plane_matrix = np.vstack((self.p2.xyz() - self.p1.xyz(),
                                  self.p3.xyz() - self.p1.xyz(),
                                  self.p4.xyz() - self.p1.xyz()))
        if round(np.linalg.det(plane_matrix).transpose(), 8) == 0:
            # Rounding here happens to prevent floating point errors from
            # interfering with the function. 8 decimals chosen arbitrarily.
            return True
        else:
            return False

    def area(self):
        """Calculate area of the face, by computing the diagonals.

           Returns: face area (float)"""
        diag13 = self.p3.xyz() - self.p1.xyz()
        diag24 = self.p4.xyz() - self.p2.xyz()
        cpdiag = np.cross(diag13, diag24)
        facearea = 0.5 * np.sqrt(np.einsum('...i,...i', cpdiag, cpdiag))
        return facearea

    def vertices(self):
        """Returns a list of the four vertices spanning the face."""
        return [self.p1, self.p2, self.p3, self.p4]

    def corners(self):
        """Synonym function of self.vertices()."""
        return self.vertices()

    def plotlist_local(self):
        """Return three lists with the x, y, and z-components of all four
            vertices in the face."""

        xlist = []
        ylist = []
        zlist = []

        for vertex in self.vertices():
            xlist.append(vertex.x)
            ylist.append(vertex.y)
            zlist.append(vertex.z)

        return xlist, ylist, zlist

    def project(self, plane='xy', new_frame=None):
        """Project a copy of the face onto a plane that is spanned by two
            axes. The projection is orthographic.

            If this face is not attached to the global coordinate frame, then
            the face returned by this function will be in terms of the local
            frame.

            TODO: A "new_frame" can be specified, which will make the
            projection elements children of this frame. Otherwise, the
            elements will be children of the global coordinate frame.

            TODO: Generalize this function for arbitrary projection plane.
            """

        # Eliminate minus sign in front of the plane:
        if plane[0] == '-':
            plane = plane[1:]

        pj_xyz = []

        # Flatten the coordinates of each vertex onto the projection plane:
        for vertex in [self.p1, self.p2, self.p3, self.p4]:
            if plane == 'xy':
                pj_x = vertex.x
                pj_y = vertex.y
                pj_z = 0
            elif plane == 'xz':
                pj_x = vertex.x
                pj_y = 0
                pj_z = vertex.z
            elif plane == 'yz':
                pj_x = 0
                pj_y = vertex.y
                pj_z = vertex.z
            else:
                raise ValueError("No valid projection plane given to "
                                 "projection method! "
                                 "Valid options: 'xy', 'xz', 'yz'")
            pj_xyz.append(Vertex([pj_x, pj_y, pj_z], parenttype="face"))

        # Create the projected face using the projected vertices.
        if not new_frame:
            pj = Face(pj_xyz[0], pj_xyz[1], pj_xyz[2], pj_xyz[3])
        else:
            pj = Face(pj_xyz[0], pj_xyz[1], pj_xyz[2], pj_xyz[3],
                      parenttype="frame")
            new_frame.add_face(pj)
        return pj

    def find_centroid(self):
        """Find the vertex centroid of the face, and then return its
           local coordinates in an numpy.ndarray.

           Returns: numpy.ndarray([xc, yc, zc])
           """
        xc = 0.25 * (self.p1.x + self.p2.x + self.p3.x + self.p4.x)
        yc = 0.25 * (self.p1.y + self.p2.y + self.p3.y + self.p4.y)
        zc = 0.25 * (self.p1.z + self.p2.z + self.p3.z + self.p4.z)
        return np.array([xc, yc, zc])

    def make_centroid(self):
        """Find the vertex centroid of the face, and turn it into a Vertex.

           Returns: Vertex([xc, yc, zc])
           """
        (xc, yc, zc) = self.find_centroid()
        return Vertex([xc, yc, zc])

    def find_perpendicular(self):
        """Find a vector perpendicular to a face (direction is ambiguous).
           TODO: Generalize the direction of the perpendicular vector.
           Currently the direction depends on how the face is defined. """

        # Find perpendicular vector to plane spanned by p12, p14.
        p12 = (self.p2.xyz() - self.p1.xyz())
        p14 = (self.p4.xyz() - self.p1.xyz())
        perpendicular = np.cross(p12, p14)

        # Then normalize it to a unit vector and return
        return perpendicular / np.linalg.norm(perpendicular)

    def readout(self, dec=4):
        """Print the vertex coordinates of all corners."""
        for vertex in [self.p1, self.p2, self.p3, self.p4]:
            print(vertex.xyz())


class Geometry:
    """In this code, a Geometry is defined as a loose collection of Faces.
       The constructor merely creates an empty Geometry object, and can
       subsequently be filled with Faces using the add_face() and
       add_faces methods.
       Some methods of this class interact with the Face objects, whilst
       others interact instead with the individual Vertices out of which
       the Faces are made, to avoid e.g. applying transformations more than
       once to the same Vertex object."""

    def __init__(self, faces=(), parenttype="global"):

        self.faces = []
        self.add_faces(faces)

        self.parenttype = "global"

        # Re-assign self.parenttype to check for validity of given parenttype
        if parenttype != "global":
            self.set_parenttype(parenttype)

    def set_parenttype(self, new_parenttype):
        """Sets own parenttype to a specified parenttype. First verifies
        specified parenttype.
        """
        if new_parenttype in ["global", "frame"]:
            self.parenttype = new_parenttype
        else:
            raise ValueError("The parenttype cannot be anything other than: \
                             'global', 'frame'!")

    def add_face(self, face: Face):
        """Add a singular face to the geometry."""
        self.faces.append(face)
        face.set_parenttype("geometry")

    def add_faces(self, faces: list):
        """Add a list of faces to the geometry."""
        for face in faces:
            if isinstance(face, Face):
                self.faces.append(face)
                face.set_parenttype("geometry")
            else:
                raise TypeError("This function can only add Face objects, \
                                but object of type {} was given!"
                                .format(type(face)))

    def vertices(self):
        vertices = []
        for face in self.faces:
            vertices.extend(face.vertices())
        return list(set(vertices))

    def translate(self, dx=0., dy=0., dz=0.):
        """Translate all the vertices in the geometry."""
        vertices = self.vertices()
        for vertex in vertices:
            vertex.translate(dx=dx, dy=dy, dz=dz)

    def rotate(self, a, b, c, cor, seq='321'):
        """Rotate all the vertices in the geometry around a point in space.
           a is Euler angle of rotation around x, etc...
               expressed in radians
           seq is a string with the rotation sequence, e.g. '321' for:
               Rz(c).Ry(b).Rx(a).vertex
           cor is the centre of rotation, which must be specified as
               a vertex, coordinate list, or coordinate numpy.ndarray
           """
        vertices = self.vertices()
        for vertex in vertices:
            vertex.rotate(a, b, c, cor=cor, seq=seq)

    def area(self):
        """Total area of all faces in the geometry.
        Note that faces have a single side, not two sides. """
        A = 0
        for face in self.faces:
            A += face.area
        return A

    def find_fcentroids(self):
        """Returns a list of the centroid coordinates of all faces currently
            in the geometry. Technically, these are vertex centroids.
            Optionally, centroids can be returned as a Numpy array by
            setting returnarray to True."""
        fcentroids = []
        for face in self.faces:
            fcentroids.append(list(face.find_centroid()))
        return fcentroids

    def find_cuboid_centroid(self):
        """Locate the centroid of a cuboid geometry.
           Cuboids have eight corners, which this method checks first.
           """
        vertices = self.vertices()
        if len(vertices) != 8:
            raise ValueError("Tried to use method find_cuboid_centroid() \
                             with a geometry containing {} vertices. \
                             All cuboids must have exactly 8 vertices. "
                             .format(len(vertices)))
        xc = 0
        yc = 0
        zc = 0
        for vertex in vertices:
            xc += vertex.x / len(vertices)
            yc += vertex.y / len(vertices)
            zc += vertex.z / len(vertices)
        return np.array([xc, yc, zc])

    def make_cuboid_centroid(self):
        """Locates the centroid of a cuboid geometry and creates a
        Vertex at this point."""
        (xc, yc, zc) = self.find_cuboid_centroid()
        return Vertex([xc, yc, zc])


class Vector:
    """Custom vector implementation. Unfinished.

    TODO: Implement translation.
    TODO: Implement rotation.
    TODO: Implement length/norm.
    TODO: Implement various numpy linalg compatibility.
    TODO: ...
    """

    def __init__(self, head: Vertex, tail: Vertex = Vertex([0, 0, 0]),
                 parenttype="global"):

        self.head = head  # Vertex
        self.tail = tail  # Vertex

        self.parenttype = "global"

        # Re-assign self.parenttype to check for validity of given parenttype
        if parenttype != "global":
            self.set_parenttype(parenttype)

        # Let child vertices know their parenttype is 'vector'
        for vertex in self.vertices():
            vertex.set_parenttype("vector")

    def set_parenttype(self, new_parenttype):
        if new_parenttype in ["global", "frame"]:
            self.parenttype = new_parenttype
        else:
            raise ValueError("The parenttype cannot be anything other than: \
                             'global', 'frame'!")

    def vertices(self):
        """Returns a list of the two vertices defining the vector."""
        return [self.head, self.tail]


class Frame:
    """Custom implentation of a 3D point expressed in cartesians.

       TODO: Clean up projection/illumination functions in a more logical
               manner.
       TODO: Add support for parenttype being "frame", in which case frames
               can be daisy-chained right up to the global frame.
    """

    def __init__(self, xyz=(0., 0., 0.), parenttype="global"):
        self.x = xyz[0]
        self.y = xyz[1]
        self.z = xyz[2]
        self.xdir = np.array([1, 0, 0])
        self.ydir = np.array([0, 1, 0])
        self.zdir = np.array([0, 0, 1])

        self.dcm = None
        self.recalculate_dcm()

        self.vertices = []
        self.faces = []
        self.geometries = []
        self.vectors = []

        self.parenttype = "global"
        # Re-assign self.parenttype to check for validity of given parenttype
        if parenttype != "global":
            self.set_parenttype(parenttype)

    def add_vertex(self, vertex: Vertex):
        self.vertices.append(vertex)
        vertex.set_parenttype = "frame"

    def remove_vertex(self, vertex: Vertex):
        self.vertices.remove(vertex)
        vertex.set_parenttype = "global"

    def add_face(self, face: Face):
        self.faces.append(face)
        face.set_parenttype = "frame"

    def remove_face(self, face: Face):
        self.faces.remove(face)
        face.set_parenttype = "global"

    def add_geometry(self, geometry):
        if isinstance(geometry, Geometry):
            self.geometries.append(geometry)
        elif isinstance(geometry, list):
            for geo in geometry:
                if isinstance(geo, Geometry):
                    self.geometries.append(geometry)
                else:
                    raise TypeError("Frame.add_geometry() argument must be a \
                                    list of geometries!")
        else:
            raise TypeError("Invalid argument provided to \
                            Frame.add_geometry()!")

    def remove_geometry(self, geometry: Geometry):
        """TODO: merge into one add_child method."""
        self.geometries.remove(geometry)
        geometry.set_parenttype = "global"

    def add_vector(self, vector: Vector):
        """TODO: merge into one add_child method."""
        self.vectors.append(vector)
        vector.set_parenttype = "frame"

    def remove_vector(self, vector: Vector):
        """TODO: merge into one add_child method."""
        self.vectors.remove(vector)
        vector.set_parenttype = "global"

    def origin(self):
        """Returns the origin of the frame.

           Returns: origin (numpy.ndarray)
        """
        return np.array([self.x, self.y, self.z])

    def recalculate_dcm(self):
        """Computes the direction cosine matrix (DCM) for conversion between
           the frame's local coordinates and the global frame.

           Returns: dcm (3-by-3 numpy.ndarray)
           """
        gx = np.array([1, 0, 0])
        gy = np.array([0, 1, 0])
        gz = np.array([0, 0, 1])
        self.dcm = np.array([
            [np.dot(gx, self.xdir),
             np.dot(gx, self.ydir),
             np.dot(gx, self.zdir)],

            [np.dot(gy, self.xdir),
             np.dot(gy, self.ydir),
             np.dot(gy, self.zdir)],

            [np.dot(gz, self.xdir),
             np.dot(gz, self.ydir),
             np.dot(gz, self.zdir)]])

    def translate(self, dx=0., dy=0., dz=0.):
        """Move the frame's origin by a relative amount."""
        self.x += dx
        self.y += dy
        self.z += dz

    def rotate(self, a, b, c, cor=None, seq='321'):
        """ Rotates the frame around a point in the global frame
        .
            a, b, c is Euler angle of rotation around local x, y, z axes
                expressed in radians

            seq is a string with the rotation sequence, e.g. '321' for:
                Rz(c).Ry(b).Rx(a).vertex
                Default rotation sequence: '321' (zyx)

            cor is the centre of rotation, which can be specified as
                a Vertex, or a list or numpy.ndarray of local coordinates.
        """

        # ==== COR argument handling ====
        # Case: COR is given as a list/ndarray of coordinates
        if isinstance(cor, list) or isinstance(cor, np.ndarray):
            if len(cor) == 3:
                cor = Vertex(list(cor))

        # Case: COR is given as a Vertex
        elif isinstance(cor, Vertex):
            pass

        # If COR is None, ignore now and handle later
        elif not cor:
            pass

        # All other cases
        else:
            raise TypeError("Error when rotating frame. Centre of rotation \
                            argument not provided correctly!")

        # ==== Define relevant rotation matrices ====
        Rx = np.array([[1, 0, 0],
                       [0, cos(a), -sin(a)],
                       [0, sin(a), cos(a)]])

        Ry = np.array([[cos(b), 0, sin(b)],
                       [0, 1, 0],
                       [-sin(b), 0, cos(b)]])

        Rz = np.array([[cos(c), -sin(c), 0],
                       [sin(c), cos(c), 0],
                       [0, 0, 1]])

        # ==== Perform rotation ====
        for axis in reversed(seq):

            # If no Centre of Rotation (cor) is given, take frame origin.
            # If not, subtract coordinates of cor first before applying
            # the rotation, and then reverse this afterwards.
            temp = (self.x, self.y, self.z)
            if cor:
                self.translate(-1 * temp[0], -1 * temp[1], -1 * temp[2])

            if axis == '1':
                self.xdir = np.dot(Rx, self.xdir)
                self.ydir = np.dot(Rx, self.ydir)
                self.zdir = np.dot(Rx, self.zdir)
            if axis == '2':
                self.xdir = np.dot(Ry, self.xdir)
                self.ydir = np.dot(Ry, self.ydir)
                self.zdir = np.dot(Ry, self.zdir)
            if axis == '3':
                self.xdir = np.dot(Rz, self.xdir)
                self.ydir = np.dot(Rz, self.ydir)
                self.zdir = np.dot(Rz, self.zdir)

            # Reverse temporary translation.
            if cor:
                self.translate(1 * temp[0], 1 * temp[1], 1 * temp[2])

        # Update the frame's direction cosine matrix:
        self.recalculate_dcm()

    def vertex_xyz_global(self, vertex):
        """Vertex coordinates in terms of parent frame."""

        # If vertex object was given:
        if isinstance(vertex, Vertex):
            xyz_local = vertex.xyz()
        elif isinstance(vertex, list) or isinstance(vertex, np.ndarray):
            if len(vertex) == 3:
                xyz_local = np.array(vertex)
            else:
                raise ValueError("This function can only convert 3D vertices,\
                                  but {} coordinates were given!"
                                 .format(len(vertex)))
        else:
            raise TypeError("Function vertex_xyz_global does not accept \
                            arguments of this type!")

        # Coordinates of local frame origin in terms of global frame
        o_local = np.array([self.x, self.y, self.z])

        # Vertex coordinates in terms of global frame
        xyz_global = np.dot(self.dcm, xyz_local) + o_local

        # return xyz_global[0], xyz_global[1], xyz_global[2]
        return xyz_global

    def geometry_vertices(self):
        geometry_vertices = []
        for geometry in self.geometries:
            geometry_vertices.extend(geometry.vertices())
        return list(set(geometry_vertices))

    def plotlist(self, face: Face):
        """Return three lists with the x, y, and z-components of all four
           vertices in the face."""

        xlist = []
        ylist = []
        zlist = []

        for vertex in face.vertices():
            xg, yg, zg = self.vertex_xyz_global(vertex)
            xlist.append(xg)
            ylist.append(yg)
            zlist.append(zg)

        return xlist, ylist, zlist

    def plotlist_xyz(self, face: Face):
        """Return a list of lists with the xyz coordinates of each vertex."""
        # If object has no parent, return plot list in terms of global frame:

        plotlist_xyz = []

        for vertex in face.vertices():
            plotlist_xyz.append(list(self.vertex_xyz_global(vertex)))

        return plotlist_xyz

    def face_perpendicular(self, face: Face, scale=1.):
        """Returns the perpendicular of a specified face in terms of the
           global coordinate frame.

           Returns: perpendicular (3-by-1 numpy.ndarray)
           """
        p1 = self.vertex_xyz_global(face.p1)
        p2 = self.vertex_xyz_global(face.p2)
        p4 = self.vertex_xyz_global(face.p4)

        perpendicular = np.cross(p2 - p1, p4 - p1)
        return perpendicular / np.linalg.norm(perpendicular)

    def geometry_perpendiculars(self, geometry: Geometry, scale=1.):
        """Returns a list with a vector perpendicular to each face in the
            geometry. If scale=1, these vectors will be unit vectors.

            TODO: Guarantee that all perpendiculars point away from the
            centre of the geometry, so that the output is not a mishmash of
            inward and outward pointing vertices. Can achieve this by taking
            inner product of:
                - vector pointing from geometry centroid to Face centroid.
                - computed perpendicular vector for given face.
            If the result is positive, the perpendicular is pointing away from
            geometry, and if result is negative, flip the perpendicular.
            Downside: Does not work well for concave geometries whose centroid
            is outside the geometry boundaries.
            """
        perpendiculars = []
        for face in geometry.faces:
            perp = self.face_perpendicular(face)
            perpendiculars.append(perp * scale)

        # Code that ensures arrows point outwards (TODO):
        # ...

        return perpendiculars

    def face_centroid(self, face: Face):
        """Find the vertex centroid of the face."""

        xc = 0
        yc = 0
        zc = 0

        for vertex in face.vertices():
            (dx, dy, dz) = self.vertex_xyz_global(vertex)
            xc += 0.25 * dx
            yc += 0.25 * dy
            zc += 0.25 * dz

        return np.array([xc, yc, zc])

    def geometry_face_centroids(self, geometry: Geometry):
        """Find the vertex centroids of all the faces in specified geometry
           in terms of the global coordinate frame.

           Returns: np.ndarray
           """
        face_centroids = []
        for face in geometry.faces:
            face_centroids.append(self.face_centroid(face))
        return np.array(face_centroids)

    def perpendiculars_plotlist(self, geometry: Geometry, scale=.1):
        """Returns a list of lists with coordinates used for plotting.
            Each item in the list is a list of six coordinates:
            - The first three coordinates indicate the xyz of the tail point.
            - The second three coordinates are point a vector from this tail
              point.

            TODO: Pass scaling as argument.
            """
        c = self.geometry_face_centroids(geometry)
        p = self.geometry_perpendiculars(geometry, scale)  # Note: 0.05

        if len(c) != len(p):
            raise ValueError("Number of centroids and perpendiculars is "
                             "somehow unequal! This should not happen.")

        # Generate plotarray, which is structured as [xyzuvw1, xyzuvw2, ...]
        plotarray = np.hstack([np.vstack(c), np.vstack(p)]).transpose()
        return plotarray.tolist()

    def illumination_vector(self, plane='xy'):
        """Generate an illumination vector based on a given plane.
           This vector essentially "simulates" where the light is coming from.
           Currently limited to directions parallel to global XYZ system.

           Returns: 3-by-1 np.ndarray
           """
        if plane == 'xy':
            iv = np.array([0, 0, 1])
        elif plane == '-xy':
            iv = np.array([0, 0, -1])
        elif plane == 'xz':
            iv = np.array([0, 1, 0])
        elif plane == '-xz':
            iv = np.array([0, -1, 0])
        elif plane == 'yz':
            iv = np.array([1, 0, 0])
        elif plane == '-yz':
            iv = np.array([-1, 0, 0])
        else:
            raise ValueError("Invalid plane given. Options: 'xy', '-xy', "
                             "'xz', '-xz', 'yz', '-yz'")
        return iv

    def illuminated_faces(self, geometry: Geometry, plane='xy'):
        """Input a geometry and a illumination plane, and this function will
            return the faces in the geometry that are illuminated.
        """
        iv = self.illumination_vector(plane=plane)

        illuminated_faces = []

        for face in geometry.faces:
            if np.dot(iv, self.face_perpendicular(face)) < 0:
                illuminated_faces.append(face)
            else:
                pass
        return illuminated_faces

    def illuminated_strength(self, face: Face, plane='xy'):
        iv = self.illumination_vector(plane=plane)
        return np.arccos(np.dot(iv, self.face_perpendicular(face))) - np.pi / 2

    def cosi_face(self, face: Face, vector):
        """Calculates the cosine of the angle between the perpendicular of a
           specified face and a specified vector.
           """
        return np.dot(vector, self.face_perpendicular(face))

    def area_projected(self, geometry: Geometry, vector):
        """Takes a geometry and a vector in a certain direction, expressed
           in the global frame. It then only looks at the faces of this
           geometry whose perpendiculars are pointing in the opposite
           direction as the specified vector. Then it computes the projected
           area of these faces with respect to the specified vector.

           In other words, if you were to look at a geometry, and the
           specified vector indicates the direction from which you observe it,
           this function returns the area of the geometry you would observe.

           returns: A_projected_total (float)
           """
        A_partial = []

        for face in geometry.faces:
            cosi = self.cosi_face(face, vector)
            if cosi < 0:
                A_partial.append(-cosi * face.area())
            else:
                A_partial.append(0)
        return sum(A_partial)

    def project_face(self, face: Face, plane='xy', new_frame=None):
        """Project a copy of the face onto a plane that is spanned by two
            axes. The projection is orthographic. Coordinates will be in terms
            of the GLOBAL frame.

            TODO: A "new_frame" can be specified, which will make the
            projection elements children of this frame. Otherwise, the
            elements will be children of the global coordinate frame.

            TODO: Generalize this function for arbitrary projection plane.
            """
        # Eliminate minus sign in front of the plane:
        if plane[0] == '-':
            plane = plane[1:]

        pj_xyz = []

        for vertex in [face.p1, face.p2, face.p3, face.p4]:

            global_coords = self.vertex_xyz_global(vertex)

            if plane == 'xy':
                pj_x = global_coords[0]
                pj_y = global_coords[1]
                pj_z = 0
            elif plane == 'xz':
                pj_x = global_coords[0]
                pj_y = 0
                pj_z = global_coords[2]
            elif plane == 'yz':
                pj_x = 0
                pj_y = global_coords[1]
                pj_z = global_coords[2]
            else:
                raise ValueError("No valid projection plane given to "
                                 "projection method! "
                                 "Valid options: 'xy', 'xz', 'yz'")
            pj_xyz.append(Vertex([pj_x, pj_y, pj_z], parenttype="face"))

        # Create the projected face using the projected vertices.
        if not new_frame:
            pj = Face(pj_xyz[0], pj_xyz[1], pj_xyz[2], pj_xyz[3],
                      parenttype="global")
        else:
            raise TypeError("new_frame functionality is not yet implemented.")
        return pj

    def illuminated_area(self, geometry: Geometry, plane='xy'):
        area = 0
        for face in self.illuminated_faces(geometry, plane):
            area += self.project_face(face, plane).area()
        return area

    def readout(self, dec=4):
        """More detailed information about Vertex.
           Use dec to define the number of decimals before rounding.
           Set dec to -1 to disable rounding.
           """
        if dec == -1:  # Set dec to -1 to disable rounding
            origin = self.origin()
        else:
            origin = np.round(self.origin(), dec)

        # children = len(self.vertices)   + \
        #            len(self.faces)      + \
        #            len(self.geometries) + \
        #            len(self.vectors)

        print("Carthesian coordinate frame.")
        print("Frame origin: {}\n".format(origin))
        # print("Children:              {}".format(children))
        # print("Associated Vertices:   {}".format(len(self.vertices)))
        # print("Associated Faces:      {}".format(len(self.faces)))
        print("Associated Geometries: {}".format(len(self.geometries)))
        # print("Associated Vectors:    {}".format(len(self.vectors)))


#%% Utilities ================================================================

def r2d(a):
    """Convert radians to degrees."""
    return a * 180 / np.pi

def d2r(a):
    """Convert degrees to radians."""
    return a * np.pi / 180

def e2q(a, b, c, seq='321'):
    """Converts Euler angles to quaternions."""
    q = np.array([
        [sin(.5 * c) * cos(.5 * b) * cos(.5 * a) - cos(.5 * c) * sin(.5 * b) * sin(.5 * a)],
        [cos(.5 * c) * sin(.5 * b) * cos(.5 * a) + sin(.5 * c) * cos(.5 * b) * sin(.5 * a)],
        [cos(.5 * c) * cos(.5 * b) * sin(.5 * a) - sin(.5 * c) * sin(.5 * b) * cos(.5 * a)],
        [cos(.5 * c) * cos(.5 * b) * cos(.5 * a) + sin(.5 * c) * sin(.5 * b) * sin(.5 * a)]
    ])
    return q

def q2e(quaternion, seq='321', returnarray=False):
    """Converts a unit quaternion to a set of euler angles."""
    q1, q2, q3, q4 = quaternion
    q1 = float(q1)
    q2 = float(q2)
    q3 = float(q3)
    q4 = float(q4)
    length = round(np.sqrt(q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4), 8)
    if length != 1:
        raise ValueError("Quaternion given is not a unit quaternion. "
                         "|q|^2 = {} != 1".format(length))

    a = np.arctan2(2 * (q4 * q3 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3))
    b = np.arcsin(2 * (q4 * q2 - q3 * q1))
    c = np.arctan2(2 * (q4 * q1 + q2 * q3), 1 - 2 * (q1 * q1 + q2 * q2))

    if returnarray is True:
        return np.array([round(a, 10), round(b, 10), round(c, 10)])
    else:
        return round(a, 10), round(b, 10), round(c, 10)

def hex2nRGB(hex_colour):
    """Converts 3, 6, or 9-length hexidecimal codes to normalized RGB values.
        The normalization is around 1, i.e. 00 -> 0.0 and FF -> 1.0

        Returns: tuple(normalized R, normalized G, normalized B)
        """
    # Strip hash off colour code:
    hex_colour = hex_colour.lstrip('#')

    # Verifying formatting (check if all hexidecimal numbers):
    if not all(c.lower() in "0123456789abcdef" for c in hex_colour):
        raise ValueError("Hex colour code '#"
                         + hex_colour
                         + "' contains illegal characters!")
    else:
        # 1-digit numbers for r/g/b:
        if len(hex_colour) == 3:
            return tuple(int(hex_colour[i:i + 1], 16) / 15 for i in (0, 1, 2))
        # 2-digit numbers for r/g/b:
        elif len(hex_colour) == 6:
            return tuple(int(hex_colour[i:i + 2], 16) / 255 for i in (0, 2, 4))
        # 3-digit numbers for r/g/b:
        elif len(hex_colour) == 9:
            return tuple(int(hex_colour[i:i + 3], 16) / 4095 for i in (0, 3, 6))
        else:
            raise ValueError("Hex colour code '#"
                             + hex_colour
                             + "' is of incorrect length!")


#%% Plotting =================================================================

def plot_vertex(axes: plt.matplotlib.axes, vertex: Vertex,
                xyz_global=None,
                vertexfill=True,  # If False, vertex will not be plotted
                vertexcolour="#000",  # Specifies the vertex colour
                vertexsize=10,  # Size of the plotted vertex
                vertexalpha=1  # Opacity of the plotted vertex
                ):
    # print("[DEBUG] Plotting {}".format(vertex))

    # Check whether vertex should be plotted:
    if vertexfill:
        # If vertex is an orphan, use its local coordinates for plotting
        # if vertex.parenttype == "global":
        if not xyz_global.all():
            x, y, z = vertex.xyz()
        # If not, check and use xyz_global coordinates provided
        else:
            # Verifying type of xyz_global
            if (isinstance(xyz_global, list) or isinstance(xyz_global, np.ndarray)):
                # Verifying length of xyz_global
                if len(xyz_global) == 3:
                    [x, y, z] = xyz_global
                else:
                    raise ValueError("Error: Tried to plot vertex using {}\
                                     coordinates, but exactly 3 are needed!\
                                     ".format(len(xyz_global)))
            # elif xyz_global == None:
            #     raise TypeError("Something went wrong with vertex plot!")
            else:
                raise TypeError("Error: Tried to plot vertex, but argument \
                                'xyz_global' was invalid type ({})!\
                                ".format(type(xyz_global)))

        axes.scatter(x, y, z,
                     c=vertexcolour, s=vertexsize, alpha=vertexalpha)


def plot_frame(axes: plt.matplotlib.axes, frame: Frame,
               # Tripod properties
               show_tripod=True,  # If False, does not plot the tripod
               tripod_scale=1,  # Sets the scale of the tripod
               plot_perpendiculars=True,  # Plot perpendiculars

               # Vertex plotting properties:
               vertexfill=True,  # If False, vertex will not be plotted
               vertexcolour="#000",  # Specifies the vertex colour
               vertexsize=10,  # Size of the plotted vertex
               vertexalpha=1,  # Opacity of the plotted vertex

               # Face plotting properties:
               linefill=True,  # If False, does not plot face lines
               linecolour="#000",  # Colour of face lines
               linewidth=2,  # Thickness of face lines
               linealpha=1,  # Opacity of face lines
               facefill=True,  # If False, does not shade the face area
               facecolour="#555",  # Colour of the face area shading
               facealpha=1,  # Opacity of the face area shading

               # Face perpendicular arrow plotting properties:
               perpfill=False,  # If True, plot perpendiculars
               perpcolour="#888",  # Specifies the perp. arrow colour
               perpscale=1,  # Size of the plotted perp. arrow
               perpalpha=0.5,  # Opacity of the plotted perp. arrow

               # Illumination:
               illumination=False,  # If True, plots illumination intensity
               ill_value=0,  # Used to plot illumination intensity
               ill_plane=None,  # If illumination is used, a plane
               #  MUST be satisfied.

               # Vector plotting properties:
               vectorfill=True,  # If False, does not plot vector arrow
               vectorcolour="#000",  # Colour of vector arrow
               vectoralpha=1,  # Opacity of vector arrow
               vectorscale=1,  # Scale the whole vector by a constant
               vectorratio=0.15  # Vector arrow length ratio
               ):
    if frame.geometries:
        """If frame contains geometries, plot these. Then find out if there
           are any leftover edges and vertices associated with the frame, plot
           these too.
           """

        for geometry in frame.geometries:

            # Plot unique vertices first:
            for vertex in geometry.vertices():
                plot_vertex(axes, vertex, frame.vertex_xyz_global(vertex),
                            vertexfill=vertexfill,
                            vertexcolour=vertexcolour,
                            vertexsize=vertexsize,
                            vertexalpha=vertexalpha)

            # Then plot faces, without plotting the vertices
            for face in geometry.faces:
                if illumination:
                    if not ill_plane:
                        raise AssertionError("If illumination is plotted, an \
                                             illumination plane MUST be \
                                             specified, e.g. \
                                             (ill_plane = 'xy')")
                    else:
                        if face not in frame.illuminated_faces(geometry,
                                                               plane=ill_plane
                                                               ):
                            ill_value = 0
                        else:
                            ill_value = frame.illuminated_strength(face,
                                                                   plane=ill_plane)

                plot_face(axes, face, frame.plotlist(face),
                          linefill=linefill,
                          linecolour=linecolour,
                          linewidth=linewidth,
                          linealpha=linealpha,
                          facefill=facefill,
                          facecolour=facecolour,
                          facealpha=facealpha,
                          vertexfill=False,
                          vertexcolour=vertexcolour,
                          vertexsize=vertexsize,
                          vertexalpha=vertexalpha,
                          illumination=illumination,
                          ill_value=ill_value)

            if perpfill:
                # First plot centroids of each plane:
                plot_geometry_perpendiculars(
                    axes,
                    frame.perpendiculars_plotlist(geometry, perpscale),
                    colour=perpcolour,
                    alpha=perpalpha
                )

    if show_tripod:
        plot_frame_tripod(axes, frame, scaling=tripod_scale)


def plot_geometry_perpendiculars(axes, perpendiculars_plotlist,
                                 colour='#888', alpha=.5, perp_scale=1.):
    """Plots the normals/perpendiculars of each face in a geometry, and
        displays them as little gray arrows.
        """
    xyzuvw = perpendiculars_plotlist

    # First plot arrow bases as dots:
    axes.scatter(xyzuvw[0], xyzuvw[1], xyzuvw[2], c=colour, s=2, alpha=alpha)

    # Then, attach plane-perpendicular arrows to these points:
    axes.quiver(xyzuvw[0][:], xyzuvw[1][:], xyzuvw[2][:],
                xyzuvw[3][:], xyzuvw[4][:], xyzuvw[5][:],
                color=colour, alpha=alpha)


def plot_face(axes: plt.matplotlib.axes, face: Face,
              plotlist: list = None,

              # Face plotting properties:
              linefill=True,  # If False, does not plot face lines
              linecolour="#000",  # Colour of face lines
              linewidth=2,  # Thickness of face lines
              linealpha=1,  # Opacity of face lines
              facefill=True,  # If False, does not shade the face area
              facecolour="#555",  # Colour of the face area shading
              facealpha=1,  # Opacity of the face area shading

              # Vertex plotting properties:
              vertexfill=True,  # If False, vertex will not be plotted
              vertexcolour="#000",  # Specifies the vertex colour
              vertexsize=10,  # Size of the plotted vertex
              vertexalpha=1,  # Opacity of the plotted vertex

              # Illumination:
              illumination=False,  # If True, plots illumination intensity
              ill_value=0  # Used to plot illumination intensity
              ):
    """ Plots an individual Face object. Check source code for kwargs.
    """

    def plotlist2plotlist_xyz(plotlist):
        """Return a list of lists with the xyz coordinates of each vertex."""
        # If object has no parent, return plot list in terms of global frame:

        p1_xyz = []
        p2_xyz = []
        p3_xyz = []
        p4_xyz = []

        for coordinate_row in plotlist:
            p1_xyz.append(coordinate_row[0])
            p2_xyz.append(coordinate_row[1])
            p3_xyz.append(coordinate_row[2])
            p4_xyz.append(coordinate_row[3])

        return [p1_xyz, p2_xyz, p3_xyz, p4_xyz]

    # Situations:
    # 1. Face is an orphan (parent == "global") and so it has to use
    #    local coordinates only
    # 2. Face is part of a geometry, which we will assume has a frame underlying
    #    it

    # print("[DEBUG] Plotting {}".format(face))

    # Case 1: Face is an orphan:
    if face.parenttype == "global":
        (xlist, ylist, zlist) = face.plotlist_local()

        # Plot individual vertices:
        if vertexfill:
            for vertex in face.vertices():
                plot_vertex(axes, vertex)

        # Plot edges that connect vertices:
        for i in range(-1, len(xlist) - 1):
            axes.plot3D([xlist[i], xlist[i + 1]],
                        [ylist[i], ylist[i + 1]],
                        [zlist[i], zlist[i + 1]],
                        linecolour, alpha=linealpha, lw=linewidth)

        # Plot the face surface:
        if facefill:
            plotlist_xyz = plotlist2plotlist_xyz(face.plotlist_local())
            fs = mp3d.art3d.Poly3DCollection([plotlist_xyz],
                                             linewidth=0)

            (nR, nG, nB) = hex2nRGB(facecolour)
            fs.set_facecolor((nR, nG, nB, facealpha))
            axes.add_collection3d(fs)

    # Case 2: Face is in a frame
    else:
        (xlist, ylist, zlist) = plotlist

        # Plot individual vertices:
        # if vertexfill:
        #     for vertex in face.vertices():
        #         plot_vertex(axes, vertex)

        # Plot edges that connect vertices:
        for i in range(-1, len(xlist) - 1):
            axes.plot3D([xlist[i], xlist[i + 1]],
                        [ylist[i], ylist[i + 1]],
                        [zlist[i], zlist[i + 1]],
                        linecolour, alpha=linealpha, lw=linewidth)

        # Plot the face surface:
        if facefill:
            plotlist_xyz = plotlist2plotlist_xyz(plotlist)
            fs = mp3d.art3d.Poly3DCollection([plotlist_xyz],
                                             linewidth=0)
            if not illumination:
                (nR, nG, nB) = hex2nRGB(facecolour)
                fs.set_facecolor((nR, nG, nB, facealpha))
                axes.add_collection3d(fs)
            else:
                value = round(0.1 + 0.9 / 1.6 * ill_value, 2)
                fs.set_facecolor((value, value, value, facealpha))
                axes.add_collection3d(fs)


def plot_illumination(axes: plt.matplotlib.axes, frame: Frame,
                      plane='xy', geometry: Geometry = None,

                      # Tripod properties
                      show_tripod=False,  # If False, does not plot the tripod
                      tripod_scale=1,  # Sets the scale of the tripod
                      plot_perpendiculars=False,  # Plot perpendiculars

                      # Vertex plotting properties:
                      vertexfill=False,  # If False, vertex will not be plotted
                      vertexcolour="FFC125",  # Specifies the vertex colour
                      vertexsize=10,  # Size of the plotted vertex
                      vertexalpha=1,  # Opacity of the plotted vertex

                      # Face plotting properties:
                      linefill=True,  # If False, does not plot face lines
                      linecolour="#FFC125",  # Colour of face lines
                      linewidth=2,  # Thickness of face lines
                      linealpha=1,  # Opacity of face lines
                      facefill=False,  # If False, does not shade the face area
                      facecolour="FFC125",  # Colour of the face area shading
                      facealpha=1,  # Opacity of the face area shading

                      # Face perpendicular arrow plotting properties:
                      perpfill=False,  # If True, plot perpendiculars
                      perpcolour="#888",  # Specifies the perp. arrow colour
                      perpscale=1,  # Size of the plotted perp. arrow
                      perpalpha=0.5,  # Opacity of the plotted perp. arrow

                      # Illumination:
                      illumination=False,  # If True, plots illumination intensity
                      ill_value=0,  # Used to plot illumination intensity

                      # Vector plotting properties:
                      vectorfill=False,  # If False, does not plot vector arrow
                      vectorcolour="#000",  # Colour of vector arrow
                      vectoralpha=1,  # Opacity of vector arrow
                      vectorscale=1,  # Scale the whole vector by a constant
                      vectorratio=0.15  # Vector arrow length ratio
                      ):
    # Gather the faces to be plotted in the geometry. If no geometry was
    #   provided (geometry=None), plot all the illuminated faces in the frame.

    if not geometry:
        geometries = frame.geometries
    else:
        geometries = list(geometry)

    illuminated_faces = []

    # Find out what faces are illuminated:
    for geo in geometries:
        faces = frame.illuminated_faces(geo, plane=plane)
        illuminated_faces.extend(faces)

    # Plot each face after projecting them in the global frame:
    for face in illuminated_faces:
        # projected_face = frame.project_face(face)
        plot_face(axes, frame.project_face(face, plane=plane),
                  linefill=linefill, linecolour=linecolour, linewidth=linewidth,
                  linealpha=linealpha, facefill=facefill, facecolour=facecolour,
                  facealpha=facealpha, vertexfill=vertexfill,
                  vertexcolour=vertexcolour, vertexsize=vertexsize,
                  vertexalpha=vertexalpha, illumination=False)

def plot_vector(axes: plt.matplotlib.axes, vector: Vector):
    """Plots an individual vector object.."""
    pass

def plot_arrow(axes, base, head, linewidth=1.5,
               alpha=1, scaling=1., alr=0.15, color="#000"):
    """Plots an rgb xyz tripod at the global origin.
        - alpha changes the opacity of the arrows. (default: 1)
        - scaling changes the length of the arrows (default: 1)
        """
    sc = scaling
    axes.quiver(base[0], base[1], base[2],
                sc * head[0], sc * head[1], sc * head[2], linewidth=linewidth,
                arrow_length_ratio=alr, color=color, alpha=alpha)

def plot_frame_tripod(axes: plt.matplotlib.axes, frame: Frame,
                      alpha=1, scaling=1.):
    """Plots an rgb xyz tripod at the global origin.
        - alpha changes the opacity of the arrows. (default: 1)
        - scaling changes the length of the arrows (default: 1)
        """
    sc = scaling
    axes.quiver(frame.x, frame.y, frame.z,
                sc * frame.xdir[0], sc * frame.xdir[1], sc * frame.xdir[2],
                arrow_length_ratio=0.15, color='red', alpha=alpha)
    axes.quiver(frame.x, frame.y, frame.z,
                sc * frame.ydir[0], sc * frame.ydir[1], sc * frame.ydir[2],
                arrow_length_ratio=0.15, color='green', alpha=alpha)
    axes.quiver(frame.x, frame.y, frame.z,
                sc * frame.zdir[0], sc * frame.zdir[1], sc * frame.zdir[2],
                arrow_length_ratio=0.15, color='blue', alpha=alpha)

def plot_global_tripod(axes: plt.matplotlib.axes, alpha=1, scaling=1.):
    """Plots an rgb xyz tripod at the global origin.
        - alpha changes the opacity of the arrows. (default: 1)
        - scaling changes the length of the arrows (default: 1)
        """
    axes.quiver(0, 0, 0, scaling, 0, 0,
                arrow_length_ratio=0.15, color='red', alpha=alpha)
    axes.quiver(0, 0, 0, 0, scaling, 0,
                arrow_length_ratio=0.15, color='green', alpha=alpha)
    axes.quiver(0, 0, 0, 0, 0, scaling,
                arrow_length_ratio=0.15, color='blue', alpha=alpha)

def plot_A_ill(A_ill: list):
    # Calculate average area
    A_avg = round(sum(A_ill) / len(A_ill), 4)

    # Set up plot
    fig_tmp = plt.figure(figsize=(10, 7))
    ax_tmp = fig_tmp.add_subplot(111)
    ax_tmp.set_ylim([0, 0.04])
    plt.title("Illuminated CubeSat area during simulation.")
    plt.xlabel("Simulation steps")
    plt.ylabel("Illuminated area [m^2]")
    plt.text(0, 0.035, "Average area: {} m^2".format(A_avg), color='r')
    # Area progression plot
    plt.plot(range(len(A_ill)), A_ill, 'k')

    # Average area plot
    plt.plot([0, len(A_ill)], [A_avg, A_avg], 'r:')