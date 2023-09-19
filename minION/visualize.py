from minION.util.IO_processor import IO_processor
import numpy as np
from bokeh.plotting import figure, show, output_notebook
from bokeh.models import ColorBar, LinearColorMapper
from bokeh.transform import transform
from bokeh.layouts import column
from bokeh.models import TabPanel, Tabs
from bokeh.layouts import layout
from bokeh.io import push_notebook

def create_heatmap_for_rbc(df):
    

    # Convert the matrix to a format suitable for Bokeh
    heatmap_matrix = heatmap_matrix[::-1]  # Reverse the matrix to align with Bokeh's coordinate system
    df = pd.DataFrame(heatmap_matrix, columns=x, index=y).stack().rename_axis(('y', 'x')).reset_index(name='value')
    
    # Create the Bokeh heatmap plot
    colors = ['#f7fbff', '#deebf7', '#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b']
    mapper = LinearColorMapper(palette=colors, low=df.value.min(), high=df.value.max())

    p = figure(title="96-Well Plate Heatmap", x_range=x, y_range=y, width=800, height=400,
               toolbar_location=None, tools="", x_axis_label="Columns", y_axis_label="Rows")
    p.rect(x="x", y="y", width=1, height=1, source=df, fill_color=transform('value', mapper))

    color_bar = ColorBar(color_mapper=mapper, location=(0, 0))
    p.add_layout(color_bar, 'right')

    layout = column(p, sizing_mode='scale_width')  # Use the Bokeh column layout for responsive width
    output_notebook()
    show(layout)

    return p, df
