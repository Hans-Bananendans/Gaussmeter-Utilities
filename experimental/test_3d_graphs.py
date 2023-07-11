"""
Simple examples demonstrating the use of GLMeshItem.
"""


import pyqtgraph as pg
import pyqtgraph.opengl as gl

app = pg.mkQApp("GLMeshItem Example")
view = gl.GLViewWidget()
view.show()
view.setWindowTitle('pyqtgraph example: GLMeshItem')
view.setCameraPosition(distance=500)


grid = gl.GLGridItem()
grid.scale(20,20,20)
view.addItem(grid)

import numpy as np

## Example 1:
## Array of vertex positions and array of vertex indexes defining faces
## Colors are specified per-face


o = [0, 0, 100]

w_x = 185
w_y = 195
w_z = 205


# coil_w = 205
# coil_h = coil_w
# coil_sp = 0.5445 * coil_w
# coil_t = 6

def make_coil(w, h=-1, t=-1, sp=-1, color=(0.5, 0.5, 0.5, 1)):
    # If h, t, sp not defined, use this as default:
    if h == -1:
        h=w
    if t == -1:
        t=w/20
    if sp == -1:
        sp=0.5*w

    verts = np.array([
        [ (w/2-t/2),  (h/2-t/2), t/2],
        [-(w/2-t/2),  (h/2-t/2), t/2],
        [-(w/2-t/2), -(h/2-t/2), t/2],
        [ (w/2-t/2), -(h/2-t/2), t/2],
        [ (w/2+t/2),  (h/2+t/2), t/2],
        [-(w/2+t/2),  (h/2+t/2), t/2],
        [-(w/2+t/2), -(h/2+t/2), t/2],
        [ (w/2+t/2), -(h/2+t/2), t/2],
        [ (w/2-t/2),  (h/2-t/2), -t/2],
        [-(w/2-t/2),  (h/2-t/2), -t/2],
        [-(w/2-t/2), -(h/2-t/2), -t/2],
        [ (w/2-t/2), -(h/2-t/2), -t/2],
        [ (w/2+t/2),  (h/2+t/2), -t/2],
        [-(w/2+t/2),  (h/2+t/2), -t/2],
        [-(w/2+t/2), -(h/2+t/2), -t/2],
        [ (w/2+t/2), -(h/2+t/2), -t/2],
    ])

    # print(verts)

    j = 4
    faces1 = []
    for i in range(j):
        # Define 4x inner triangles
        faces1.append([i, (i+1)%j, i+j])
        # Define 4x outer triangles
        faces1.append([i, i+j, j+(i-1)%j])

    # Define 8x triangles with z-offset
    faces2 = [[coor + 2*j for coor in vertex] for vertex in faces1]

    # Define connecting side triangles:
    faces3 = []
    for i in range(j):
        # Define 4x upper side triangles
        faces1.append([i+j, j+(i+1)%j, i+3*j])
        # Define 4x lower side triangles
        faces1.append([i+j, i+3*j, 3*j+(i-1)%j])

    # Concatenate and cast into Numpy array
    faces = np.array([faces1+faces2+faces3])

    # print(faces)

    # Make all faces the same colour
    colors = np.array([[color[0], color[1], color[2], color[3]]]*len(faces))

    # print(colors)

    return verts, faces, colors


def make_coil2(w, h=-1, t=-1, sp=-1, color=(0.5, 0.5, 0.5, 1), offset=(0, 0, 0)):
    # If h, t, sp not defined, use this as default:
    if h == -1:
        h=w
    if t == -1:
        t=w/60
    if sp == -1:
        sp=0.5*w

    verts = np.array([
        [ (w/2-t/2),  (h/2-t/2), t/2],
        [-(w/2-t/2),  (h/2-t/2), t/2],
        [-(w/2-t/2), -(h/2-t/2), t/2],
        [ (w/2-t/2), -(h/2-t/2), t/2],
        [ (w/2+t/2),  (h/2+t/2), t/2],
        [-(w/2+t/2),  (h/2+t/2), t/2],
        [-(w/2+t/2), -(h/2+t/2), t/2],
        [ (w/2+t/2), -(h/2+t/2), t/2],
        [ (w/2-t/2),  (h/2-t/2), -t/2],
        [-(w/2-t/2),  (h/2-t/2), -t/2],
        [-(w/2-t/2), -(h/2-t/2), -t/2],
        [ (w/2-t/2), -(h/2-t/2), -t/2],
        [ (w/2+t/2),  (h/2+t/2), -t/2],
        [-(w/2+t/2),  (h/2+t/2), -t/2],
        [-(w/2+t/2), -(h/2+t/2), -t/2],
        [ (w/2+t/2), -(h/2+t/2), -t/2],
    ])


    # print(verts)


    faces1 = [
        [0,5,1],
        [0,4,5],
        [1,6,2],
        [1,5,6],
        [2,7,3],
        [2,6,7],
        [3,4,0],
        [3,7,4],
    ]

    # # The following flips the normals the wrong way, so do not use.
    # j=4
    # faces2 = [[coor + 2*j for coor in vertex] for vertex in faces1]

    faces2 = [
        [8,13,12],
        [8,9,13],
        [9,14,13],
        [9,10,14],
        [10,15,14],
        [10,11,15],
        [11,12,15],
        [11,8,12],
    ]


    faces3 = [
        [7,6,14],
        [7,14,15],

        [4,7,15],
        [4,15,12],

        [5,4,12],
        [5,12,13],

        [6,5,14],
        [5,14,13],
    ]

    faces4 = [
        [11,10,2],
        [11,2,3],

        [8,11,3],
        [8,3,0],

        [9,8,0],
        [9,0,1],

        [10,9,1],
        [10,1,2],
    ]

    faces = np.array(faces1+faces2+faces3+faces4)
    # print(faces)

    # print(faces)
    # Make all faces the same colour
    colors = np.array([[color[0], color[1], color[2], color[3]]]*len(faces))

    # print(colors)

    mesh = gl.GLMeshItem(vertexes=verts, faces=faces, faceColors=colors)
    mesh.translate(offset[0], offset[1], offset[2])
    return mesh


# verts, faces, colors = make_coil2(coil_w, coil_h, coil_t, coil_sp)
# verts, faces, colors = make_coil2(coil_w)

metal_white = (1.0, 1.0, 1.0, 1)
metal_gray = (0.5, 0.5, 0.5, 1)

metal_red = (1.0, 0.85, 0.85, 1)
metal_green = (0.85, 1.0, 0.85, 1)
metal_blue = (0.85, 0.85, 1.0, 1)

# Mesh item will automatically compute face normals.
mesh_x_p = make_coil2(w_x, offset=( 0.5445*w_x/2,0,o[2]), color=metal_red)
mesh_x_n = make_coil2(w_x, offset=(-0.5445*w_x/2,0,o[2]), color=metal_red)
mesh_y_p = make_coil2(w_y, offset=(0, 0.5445*w_y/2,o[2]), color=metal_green)
mesh_y_n = make_coil2(w_y, offset=(0,-0.5445*w_y/2,o[2]), color=metal_green)
mesh_z_p = make_coil2(w_z, offset=(0,0, 0.5445*w_z/2+o[2]), color=metal_blue)
mesh_z_n = make_coil2(w_z, offset=(0,0,-0.5445*w_z/2+o[2]), color=metal_blue)


mesh_x_p.rotate(90, 0,1,0, local=True)
mesh_x_n.rotate(90, 0,1,0, local=True)
mesh_y_p.rotate(90, 1,0,0, local=True)
mesh_y_n.rotate(90, 1,0,0, local=True)

# mesh1 = gl.GLMeshItem(vertexes=verts, faces=faces, faceColors=colors,
#                       smooth=True, drawEdges=True)
# mesh1.translate(coil_offset[0], coil_offset[1], coil_offset[2])

for mesh in (mesh_x_p, mesh_x_n, mesh_y_p, mesh_y_n, mesh_z_p, mesh_z_n):
    mesh.opts["drawEdges"] = True
    mesh.opts["smooth"] = True
    mesh.setGLOptions('opaque')
    mesh.setShader('shaded')
    view.addItem(mesh)

# mesh1.setGLOptions('opaque')
# mesh1.setShader('viewNormalColor')
# view.addItem(mesh1)

# print(mesh1._facenormals)

if __name__ == '__main__':
    pg.exec()
