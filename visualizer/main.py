import pandas as pd

from bokeh.layouts import column, row
from bokeh.models import Select
from bokeh.palettes import Spectral5, Viridis5
from bokeh.plotting import curdoc, figure
from bokeh.sampledata.autompg import autompg_clean as m

# constants
num_ch4_a3 = 2.69015E-05 # from methane-comparison.xlsx
epsilon_max = 500

# read data and cleanup
m = pd.read_csv('pm2.csv')
m['radius_plot'] = (1/3)*m['volume']**(1/2)
m['ch4_uc'] = m.absolute_volumetric_loading  * (num_ch4_a3 * m.volume)
del m['b']
del m['c']
del m['Unnamed: 0']


SIZES = list(range(6, 22, 3))
COLORS = Spectral5
N_SIZES = len(SIZES)
N_COLORS = len(COLORS)

columns = sorted(m.columns)

def create_figure():
    xs = m[x.value].values
    ys = m[y.value].values
    x_title = x.value.title()
    y_title = y.value.title()

    kw = dict()
    kw['title'] = "%s vs %s" % (x_title, y_title)

    p = figure(plot_height=600, plot_width=800, tools='pan,box_zoom,hover,reset', **kw)
    p.xaxis.axis_label = x_title
    p.yaxis.axis_label = y_title

    sz = 9
    if size.value != 'None':
        if len(set(m[size.value])) > N_SIZES:
            groups = pd.qcut(m[size.value].values, N_SIZES, duplicates='drop')
        else:
            groups = pd.Categorical(m[size.value])
        sz = [SIZES[xx] for xx in groups.codes]

    c = "#31AADE"
    if color.value != 'None':
        if len(set(m[color.value])) > N_COLORS:
            groups = pd.qcut(m[color.value].values, N_COLORS, duplicates='drop')
        else:
            groups = pd.Categorical(m[color.value])
        c = [COLORS[xx] for xx in groups.codes]

    p.circle(x=xs, y=ys, color=c, size=sz, line_color="white", alpha=0.6, hover_color='white', hover_alpha=0.5)

    return p


def update(attr, old, new):
    layout.children[1] = create_figure()


x = Select(title='X-Axis', value='void_fraction_geo', options=columns)
x.on_change('value', update)

y = Select(title='Y-Axis', value='absolute_volumetric_loading', options=columns)
y.on_change('value', update)

size = Select(title='Size', value='None', options=['None'] + columns)
size.on_change('value', update)

color = Select(title='Color', value='None', options=['None'] + columns)
color.on_change('value', update)

controls = column(x, y, color, size, width=200)
layout = row(controls, create_figure())

curdoc().add_root(layout)
curdoc().title = "Crossfilter"
