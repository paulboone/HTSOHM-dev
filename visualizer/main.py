import numpy as np
import pandas as pd

from bokeh.layouts import column, row
from bokeh.models import Select, ColumnDataSource, ColorBar
from bokeh.palettes import Spectral5, Viridis5, Viridis8
from bokeh.plotting import curdoc, figure
from bokeh.sampledata.autompg import autompg_clean as m
from bokeh.transform import linear_cmap

from pseudomaterial_render import show_pseudomaterial

# constants
num_ch4_a3 = 2.69015E-05 # from methane-comparison.xlsx
epsilon_max = 500

# read data and cleanup
m = pd.read_csv('pm2.csv')
m['volume_plotsize'] = (1/3)*m['volume']**(1/2)
m['atom_sites_plotsize'] = m.atom_sites * 4
m['ch4_uc'] = m.absolute_volumetric_loading  * (num_ch4_a3 * m.volume)
m['epsilon_density_log'] = np.log(m.epsilon_density)
del m['b']
del m['c']
del m['Unnamed: 0']
m_source = ColumnDataSource(m)

atoms = pd.read_csv('atoms.csv')

columns = sorted(m.columns)
columns.remove('volume_plotsize')
columns.remove('atom_sites_plotsize')

print(m.head())

colormap_overrides = {
    'atom_sites': dict(palette=Viridis8),
    'site_distribution': dict(palette=Viridis5, low=0, high=0.5)
    # 'epsilon_density': dict(palette=Viridis5, low=0, high=0.5)
}

def create_figure():
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

    p = figure(plot_height=1000, plot_width=1000,
                tooltips=tooltips, tools=["tap", "hover", "box_select", "reset"],
                title=("%s vs %s" % (y.value, x.value)))
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

def update(attr, old, new):
    layout.children[1] = create_figure()
    print('layout updated')


x = Select(title='X-Axis', value='void_fraction_geo', options=columns)
x.on_change('value', update)

y = Select(title='Y-Axis', value='absolute_volumetric_loading', options=columns)
y.on_change('value', update)

size = Select(title='Size', value='atom_sites', options=['None'] + columns)
size.on_change('value', update)

color = Select(title='Color', value='site_distribution', options=['None'] + columns)
color.on_change('value', update)

controls = column(x, y, color, size, width=200)
layout = row(controls, create_figure())

curdoc().add_root(layout)
curdoc().title = "Pseudomaterial Visualizer"
