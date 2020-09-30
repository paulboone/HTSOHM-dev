from matplotlib import cm
from matplotlib.colors import rgb2hex
from pythreejs import Geometry, Line, LineBasicMaterial, SphereGeometry, MeshLambertMaterial, \
                      Mesh, PerspectiveCamera, DirectionalLight, AmbientLight, Scene, \
                      OrbitControls, Renderer
import pandas as pd


def pseudomaterial_render(atoms):
    c = cm.get_cmap("plasma")


    scale_axis_vertices = [[-1,-1,-1], [9, -1, -1]]
    scale_line_geom = Geometry(vertices=scale_axis_vertices, colors = ['black']*len(scale_axis_vertices))
    scale_lines = Line(geometry=scale_line_geom,
                 material=LineBasicMaterial(linewidth=50, vertexColors='VertexColors'),
                 type='LinePieces',
                )


    a = atoms.a[1]
    cube_vertices = [[0, 0, 0],
                     [a, 0, 0], [a, 0, a], [a, 0, 0],
                     [a, a, 0], [a, a, a], [a, a, 0],
                     [0, a, 0], [0, a, a], [0, a, 0],
                     [0, 0, 0], [0, 0, a], [a, 0, a], [a, a, a], [0, a, a], [0, 0, a]]

    linesgeom = Geometry(vertices=cube_vertices, colors = ['black']*len(cube_vertices))

    lines = Line(geometry=linesgeom,
                 material=LineBasicMaterial(linewidth=5, vertexColors='VertexColors'),
                 type='LinePieces',
                )

    balls = []
    for p in atoms.itertuples():
        positions = (p.x * p.a, p.y * p.a, p.z * p.a)
        new_ball = Mesh(geometry=SphereGeometry(radius=p.sigma, widthSegments=32, heightSegments=24),
                    material=MeshLambertMaterial(color=rgb2hex(c(p.epsilon_norm))),
                    position=positions)
        balls.append(new_ball)

#             [scale*2*a,scale*2*a,scale*2*a]
    camera = PerspectiveCamera(position=[25,a,a], up=[0, 1, 0],
                children=[DirectionalLight(color='white', position=[3, 5, 1], intensity=0.5)])

    scene = Scene(children=[scale_lines, lines, *balls, camera, AmbientLight(color='#777')])

    renderer = Renderer(camera=camera,
                        scene=scene,
                        controls=[OrbitControls(controlling=camera,
                                                target=[a,a,a])])


    return renderer


def scale_2x2(df):
    scale = [(0,0,0), (1,0,0), (0,1,0), (0,0,1), (1,1,0), (1,0,1), (0,1,1), (1,1,1)]

    p2 = pd.DataFrame()
    for offset in scale:
        poffset = df.copy()
        poffset[['x','y','z']] += offset
        p2 = pd.concat([p2,poffset], ignore_index=True)

    return p2

def show_pseudomaterial(id, atoms, epsilon_max):
    p = atoms[atoms.material_id==id]
    p = p.assign(epsilon_norm=p.epsilon / epsilon_max)
    scaled_p = scale_2x2(p)
    return pseudomaterial_render(scaled_p)
