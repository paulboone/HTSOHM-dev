import pandas as pd

from bokeh.layouts import column, row
from bokeh.models import Select, ColumnDataSource, ColorBar
from bokeh.palettes import Spectral5, Viridis5, Viridis8
from bokeh.plotting import curdoc, figure
from bokeh.sampledata.autompg import autompg_clean as m
from bokeh.transform import linear_cmap

# constants
num_ch4_a3 = 2.69015E-05 # from methane-comparison.xlsx
epsilon_max = 500

# read data and cleanup
m = pd.read_csv('pm2.csv')
m['volume_plotsize'] = (1/3)*m['volume']**(1/2)
m['ch4_uc'] = m.absolute_volumetric_loading  * (num_ch4_a3 * m.volume)
del m['b']
del m['c']
del m['Unnamed: 0']

m_source = ColumnDataSource(m)

SIZES = list(range(6, 22, 3))
COLORS = Spectral5
N_SIZES = len(SIZES)
N_COLORS = len(COLORS)

columns = sorted(m.columns)
columns.remove('volume_plotsize')

print(m.head())

TOOLTIPS = [
    ("id", "@id"),
    ("CH4 / UC", "@ch4_uc"),
    ("lattice", "@a")
]

colormap_overrides = {
    'atom_sites': dict(palette=Viridis8, low=1, high=8),
    'site_distribution': dict(palette=Viridis5, low=0, high=0.5)
}

def create_figure():
    print("creating figure with x = %s, y = %s, color = %s, size = %s" % (x.value, y.value, color.value, size.value))
    x_title = x.value.title()
    y_title = y.value.title()

    kw = dict()
    kw['title'] = "%s vs %s" % (y_title, x_title)

    p = figure(plot_height=1000, plot_width=1000,
                tooltips=TOOLTIPS, tools=["tap", "hover", "box_select", "reset"],
                **kw)
    p.xaxis.axis_label = x_title
    p.yaxis.axis_label = y_title

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
            colormap_args = dict(palette=Viridis5, low=m[color.value].min(), high=m[color.value].max())
        mapper = linear_cmap(field_name=color.value, **colormap_args)
        c = mapper

    p.circle(x=x.value, y=y.value, color=c, size=sz, line_color=c, alpha=0.4,
            hover_color='white', hover_alpha=0.7,
            source=m_source)

    if mapper:
        color_bar = ColorBar(color_mapper=mapper['transform'], width=8,  location=(0,0))
        p.add_layout(color_bar, 'right')

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
