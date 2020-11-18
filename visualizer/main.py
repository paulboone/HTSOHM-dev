from glob import glob
from os.path import basename, splitext

import numpy as np
import pandas as pd

from bokeh.layouts import column, row
from bokeh.models import Select, ColumnDataSource, ColorBar
from bokeh.palettes import Viridis5, Viridis8
from bokeh.plotting import curdoc, figure
from bokeh.transform import linear_cmap
from bokeh.models.widgets import Slider

# constants
num_ch4_a3 = 2.69015E-05 # from methane-comparison.xlsx

# read data and cleanup
def load_data(path):
    print("loading new data from %s" % path)
    m = pd.read_csv(path)
    m.rename(columns={'a'                                   : 'lattice size',
                      'atom_sites'                          : 'num atoms',
                      'bin12'                               : 'bin: void fraction',
                      'bin13'                               : 'bin: methane loading',
                      'number_density'                      : "number density",
                      'total_epsilon'                       : "total epsilon",
                      'epsilon_density'                     : "epsilon density",
                      'void_fraction'                       : "void fraction raspa",
                      'void_fraction_geo'                   : "void fraction geo",
                      'absolute_volumetric_loading'         : "absolute volumetric loading",
                      'absolute_volumetric_loading_error'   : "absolute volumetric loading error",
                      'max_pair_distance'                   : "max pair distance"
                     }, inplace=True)

    m['volume_plotsize'] = (1/3)*m['volume']**(1/2)
    m['num atoms_plotsize'] = m["num atoms"] * 4
    m['CH4 / uc'] = m["absolute volumetric loading"]  * (num_ch4_a3 * m.volume)
    m['epsilon density [log]'] = np.log(m["epsilon density"])
    m['number density [log]'] = np.log(m["number density"])

    del m['b']
    del m['c']

    m_source = ColumnDataSource(m)

    # atoms = pd.read_csv('atoms.csv')

    columns = sorted(m.columns)
    columns.remove('volume_plotsize')
    columns.remove('num atoms_plotsize')
    columns.remove('id')
    columns.remove('parent_id')

    print(m.head())
    print(columns)
    return m, m_source, columns

colormap_overrides = {
    'num atoms': dict(palette=Viridis8),
    'max pair distance': dict(palette=Viridis5, low=0, high=0.5)
    # 'epsilon_density': dict(palette=Viridis5, low=0, high=0.5)
}

range_defaults = {
    "absolute volumetric loading": (0,800),
    "void fraction geo": (0,1)
}

def create_figure(m, m_source, columns):
    print("creating figure with x = %s, y = %s, color = %s, size = %s" % (x.value, y.value, color.value, size.value))

    tooltips = [
        ("id", "@id"),
        (x.value, "@" + x.value),
        (y.value, "@" + y.value)
    ]
    if color.value != 'None':
        tooltips += [(color.value, "@" + color.value)]
    if size.value != 'None':
        tooltips += [(size.value, "@" + size.value)]

    x_range = range_defaults[x.value] if x.value in range_defaults else None
    y_range = range_defaults[y.value] if y.value in range_defaults else None

    p = figure(plot_height=800, plot_width=800, x_range=x_range, y_range=y_range,
                tooltips=tooltips, tools=["tap", "hover", "box_select", "reset", "save"],
                title=("%s: %s vs %s" % (data.value, y.value, x.value)))
    p.xaxis.axis_label = x.value
    p.yaxis.axis_label = y.value

    sz = 8
    print("size.value = '%s'" % size.value)
    if size.value != 'None':
        if (size.value + "_plotsize") in m:
            sz = size.value + "_plotsize"
        else:
            sz = size.value
        print(sz)

    mapper = None
    c = "#31AADE"
    if color.value != 'None':
        if color.value in colormap_overrides:
            colormap_args = colormap_overrides[color.value]
        else:
            colormap_args = dict(palette=Viridis5)

        if 'low' not in colormap_args:
            colormap_args['low'] = m[color.value].min()
        if 'high' not in colormap_args:
            colormap_args['high'] = m[color.value].max()

        print(color.value, colormap_args)

        mapper = linear_cmap(field_name=color.value, **colormap_args)
        c = mapper

    p.circle(x=x.value, y=y.value, color=c, size=sz, line_color=c, alpha=0.4,
            hover_color='white', hover_alpha=0.7,
            source=m_source)

    if mapper:
        color_bar = ColorBar(color_mapper=mapper['transform'], width=8,  location=(0,0))
        p.add_layout(color_bar, 'right')

    return p

def slider_on_change(attr, old, gen):
    m2 = m[m.generation <= gen]
    m2_source = ColumnDataSource(m2)
    layout.children[1] = create_figure(m2, m2_source, columns)
    print('generation updated')

def update_dataset(attr, old, new):
    print("loading dataset: ", new)
    data.options = data_files[new]
    data.value = default_data_file[new]
    print("dataset updated")

def update_data(attr, old, new):
    global m, m_source, columns
    m, m_source, columns = load_data("./data/%s/%s.csv" % (dataset.value, new))
    layout.children[1] = create_figure(m, m_source, columns)
    print('data and layout updated')

def update(attr, old, new):
    layout.children[1] = create_figure(m, m_source, columns)
    print('layout updated')


datasets = ["parameter-explorations", "degrees-of-freedom", "random-long"]
data_files = {k:sorted([splitext(basename(f))[0] for f in glob("./data/%s/*.csv" % k)]) for k in datasets}
default_data_file = {
    "parameter-explorations": "reference baseline",
    "degrees-of-freedom": "IME   1site",
    "random-long": "random 500K materials"
}

m, m_source, columns = load_data("./data/parameter-explorations/reference baseline.csv")

dataset = Select(title='Dataset', value="parameter-explorations", options=datasets)
dataset.on_change('value', update_dataset)

data = Select(title='Data source', value=default_data_file['parameter-explorations'], options=data_files['parameter-explorations'])
data.on_change('value', update_data)

x = Select(title='X-Axis', value='void fraction geo', options=columns)
x.on_change('value', update)

y = Select(title='Y-Axis', value='absolute volumetric loading', options=columns)
y.on_change('value', update)

size = Select(title='Size', value='lattice size', options=['None'] + columns)
size.on_change('value', update)

color = Select(title='Color', value='num atoms', options=['None'] + columns)
color.on_change('value', update)

slider = Slider(start=0, end=500, value=500, step=50, title="Generation")
slider.on_change('value', slider_on_change)

controls = column(dataset, data, x, y, color, size, slider, width=200)
layout = row(controls, create_figure(m, m_source, columns))


curdoc().add_root(layout)
curdoc().title = "Pseudomaterial Visualizer"
