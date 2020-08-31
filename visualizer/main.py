from glob import glob

import numpy as np
import pandas as pd

from bokeh.layouts import column, row
from bokeh.models import Select, ColumnDataSource, ColorBar
from bokeh.palettes import Spectral5, Viridis5, Viridis8
from bokeh.plotting import curdoc, figure
from bokeh.sampledata.autompg import autompg_clean as m
from bokeh.transform import linear_cmap
from bokeh.models.widgets import Slider

from pseudomaterial_render import show_pseudomaterial

# constants
num_ch4_a3 = 2.69015E-05 # from methane-comparison.xlsx
epsilon_max = 500
data_files = sorted(glob("./data/*.csv"))

# read data and cleanup
def load_data(path):
    print("loading new data from %s" % path)
    m = pd.read_csv(path)
    m.rename(columns={'void_fraction': 'void_fraction_raspa'}, inplace=True)

    m['volume_plotsize'] = (1/3)*m['volume']**(1/2)
    m['atom_sites_plotsize'] = m.atom_sites * 4
    m['ch4_uc'] = m.absolute_volumetric_loading  * (num_ch4_a3 * m.volume)
    m['epsilon_density_log'] = np.log(m.epsilon_density)
    m['number_density_log'] = np.log(m.number_density)

    del m['b']
    del m['c']
    # del m['Unnamed: 0']

    m_source = ColumnDataSource(m)

    # atoms = pd.read_csv('atoms.csv')

    columns = sorted(m.columns)
    columns.remove('volume_plotsize')
    columns.remove('atom_sites_plotsize')
    columns.remove('id')
    columns.remove('parent_id')

    print(m.head())
    print(columns)
    return m, m_source, columns


m, m_source, columns = load_data("./data/reference.csv")

colormap_overrides = {
    'atom_sites': dict(palette=Viridis8),
    'max_pair_distance': dict(palette=Viridis5, low=0, high=0.5)
    # 'epsilon_density': dict(palette=Viridis5, low=0, high=0.5)
}

range_defaults = {
    "absolute_volumetric_loading": (0,800),
    "void_fraction_geo": (0,1)
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

    #
    # def plot_on_selection(attr, old, new):
    #     print(attr, old, new)
    #     # for i in new:
    #     #     show_pseudomaterial(m.iloc[i].id)
    #     #
    #     if len(new) > 0:
    #         renderer = show_pseudomaterial(m.iloc[0].id, atoms, epsilon_max)
    #         print(renderer)

    # m_source.selected.on_change('indices', plot_on_selection)


    return p

def slider_on_change(attr, old, gen):
    m2 = m[m.generation <= gen]
    m2_source = ColumnDataSource(m2)
    layout.children[1] = create_figure(m2, m2_source, columns)
    print('generation updated')


def update_data(attr, old, new):
    global m, m_source, columns

    m, m_source, columns = load_data(data.value)
    layout.children[1] = create_figure(m, m_source, columns)
    print('data and layout updated')


def update(attr, old, new):
    layout.children[1] = create_figure(m, m_source, columns)
    print('layout updated')

data = Select(title='Data source', value='./data/reference.csv', options=data_files)
data.on_change('value', update_data)

x = Select(title='X-Axis', value='void_fraction_geo', options=columns)
x.on_change('value', update)

y = Select(title='Y-Axis', value='absolute_volumetric_loading', options=columns)
y.on_change('value', update)

size = Select(title='Size', value='volume', options=['None'] + columns)
size.on_change('value', update)

color = Select(title='Color', value='atom_sites', options=['None'] + columns)
color.on_change('value', update)

slider = Slider(start=0, end=500, value=500, step=50, title="Generation")
slider.on_change('value', slider_on_change)

controls = column(data, x, y, color, size, slider, width=200)
layout = row(controls, create_figure(m, m_source, columns))


curdoc().add_root(layout)
curdoc().title = "Pseudomaterial Visualizer"
