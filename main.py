import json

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.express as px
from fermi_surface import fermi_surface as fs

fs1 = fs('example/Al_25.in', vasprun='example/Al.xml', fermi=8.34)
fs1.patch_triangle()
fs1.move_triangle_bz()
fig = fs1.draw()
fig2 = fs1.draw_hist()
#fig = px.line(
#    x=["a","b","c"], y=[1,3,2], 
#    title="sample figure", height=325
#)

external_stylesheets=['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
        html.Div(className='row', children=[
        html.Div([dcc.Graph(id="graph", figure=fig),], className='four columns'),
        html.Div([dcc.Graph(id="graph2", figure=fig2),], className='four columns')]),
        html.Div([dcc.Slider(id='my-slider', min=8.0, max=9.1, step=0.1, value=8.34,),
                  html.Div(id='slider-output-container')])
])

@app.callback(
    Output("graph", "figure"), 
    [Input('my-slider', "value")])
def display_structure(value):
    return fs1.draw()

@app.callback(
    dash.dependencies.Output('slider-output-container', 'children'),
    [dash.dependencies.Input('my-slider', 'value')])
def update_output(value):
    return 'fermi energy is "{}"'.format(value)



app.run_server(debug=False, port=8080, host='0.0.0.0')
