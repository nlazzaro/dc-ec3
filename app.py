import json

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from dash.exceptions import PreventUpdate

import plotly.graph_objs as go

import pandas as pd
import numpy as np

from scipy.sparse import csr_matrix

# when a gene is selected from the graph or textbox,
# both graphs should color its cnt_normalized and the textbox should display the name;
# when the gene name is removed from the textbox or a click happens in a white space
# of the genes graph, everything should be back to default view and the textbox should be cleared


# for huge dropbox-lists: https://community.plot.ly/t/dynamic-options-for-drop-down/7550

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

strain_time = pd.read_csv('data/fc_strainVStime.csv')
strain_time['organism'] = [ 'mouse' if gene[:3] == 'ENS' else 'toxo' for gene in strain_time['Ensembl.Gene.ID'] ]

# some 'Associated.Gene.Name' in the strain vs time table are set to null;
# we set those to the corrisponding 'Ensembl.Gene.ID' value
strain_time.loc[
    strain_time['Associated.Gene.Name'].isnull(), 'Associated.Gene.Name'
    ] = strain_time[ strain_time['Associated.Gene.Name'].isnull() ]['Ensembl.Gene.ID']

glm = pd.read_csv('data/glmproj.csv')
glm.columns = ['glmproj_cell','x', 'y']

genes_acronyms = pd.read_csv('data/ensembl2genename.csv')
genes_acronyms = { v[0]: v[1] for v in genes_acronyms[ genes_acronyms.columns ].values }

# reads = pd.read_csv('data/cnt_normalized.csv')

cnt_normalized = {}
i=0
with open('data/cnt_normalized.csv', 'r') as file:
    for line in file:
        line = line.strip().split(',')
        i+=1
        if i>1:
            cnt_normalized[line[0].replace('"', '')] = csr_matrix(np.array(line[1:]).astype(float))
        else:
            cnt_columns = [ cell.replace('"', '') for cell in line[1:]]

tsne = pd.read_csv('data/tsne_coordinates.csv')

condition = pd.read_csv('data/condition.csv')
condition = pd.merge(condition, tsne, left_on='cellname', right_on=tsne.columns[0])
condition = condition[ condition['strain'] != 'ignore' ]

condition = pd.merge(condition, glm, left_on='cellname', right_on='glmproj_cell')

condition['shape'] = 'circle'
condition.loc[condition['strain'] == 'LDM', 'shape'] = 'triangle-up'
condition.loc[condition['strain'] == 'PTG', 'shape'] = 'triangle-down'

condition.loc[condition['strain'] == 'LDM', 'shape'] = 'diamond'
condition.loc[condition['strain'] == 'PTG', 'shape'] = 'square'

genes_dropdown_options = { v[0]: v[1] for v in strain_time[ ['Ensembl.Gene.ID', 'Associated.Gene.Name'] ].values }

'''
import resource
print('max ram consumption', int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000))
'''

for gene in cnt_normalized:
    if gene not in genes_dropdown_options:
        genes_dropdown_options[gene] = gene

genes_options = list(genes_dropdown_options)

for gene in list(genes_dropdown_options):
    if genes_dropdown_options[gene][:4] in ['ENSM', 'TGME']:
        if gene in genes_acronyms:
            genes_dropdown_options[gene] = genes_acronyms[gene]
        else:
            del genes_dropdown_options[gene]


genes_scatter_plot = go.Figure(
        data=[],
        layout=go.Layout(
            title='Genes Projection',
            showlegend=True,
            legend=go.layout.Legend( x=0, y=1.0 ),
            margin=go.layout.Margin(l=40, r=0, t=40, b=30),
            xaxis={ 'title': 'FC strain' },
            yaxis={ 'title': 'FC time' }
        )
    )

for organism in strain_time['organism'].unique():

    organism_df = strain_time[ strain_time['organism'] == organism ]

    genes_scatter_plot.add_trace(

        go.Scatter(
            mode='markers',
            x=organism_df['fc_strain'],
            y=organism_df['fc_time'],
            text=organism_df['Associated.Gene.Name'],
            opacity=0.8,
            marker=dict(
                size=10,
                line={'width': 0.5, 'color': 'white'},
                # color = strain_time['organism'],
                # color='Blue',
                # symbol=strain_df['shape'].values[0]
            ),
            name=organism
        )
    )

cells_scatter_plot = go.Figure(
        data=[],
        layout=go.Layout(
            # title='US Export of Plastic Scrap',
            showlegend=True,
            legend=go.layout.Legend( x=1, y=0 ),
            margin=go.layout.Margin(l=40, r=0, t=40, b=30),
            xaxis={ 'title': 'strain' },
            yaxis={ 'title': 'time' }
        )
    )

for strain in condition['strain'].unique():

    strain_df = condition[ condition['strain'] == strain ]

    # Add first scatter trace with medium sized markers
    cells_scatter_plot.add_trace(
        go.Scatter(
            mode='markers',
            x=strain_df['x'],
            y=strain_df['y'],
            text= [ v[0] + ', time: ' + str(v[1]) for v in strain_df[['cellname', 'time']].values ],
            opacity=0.7,
            marker=dict(
                # color='Blue',
                size=10,
                line={'width': 0.5, 'color': 'white'},
                symbol=strain_df['shape'].values[0]
            ),
            name=strain
        )
    )


cells_scatter_tsne = go.Figure(
    data=[],
    layout=go.Layout(
        # title='US Export of Plastic Scrap',
        showlegend=True,
        legend=go.layout.Legend(
            x=1.0,
            y=0
        ),
        margin=go.layout.Margin(l=40, r=0, t=40, b=30),
        xaxis={ 'title': 'tSNE 1' },
        yaxis={ 'title': 'tSNE2' }
    )
)

for strain in condition['strain'].unique():

    strain_df = condition[ condition['strain'] == strain ]

    # Add first scatter trace with medium sized markers
    cells_scatter_tsne.add_trace(
        go.Scatter(
            mode='markers',
            x=strain_df['tSNE_1'],
            y=strain_df['tSNE_2'],
            text= [ v[0] + ', time: ' + str(v[1]) for v in strain_df[['cellname', 'time']].values ],
            opacity=0.7,
            marker=dict(
                # color='Blue',
                size=10,
                line= {'width': 0.5, 'color': 'white'},
                symbol=strain_df['shape'].values[0]
            ),
            name=strain
        )
    )


app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}

app.layout = html.Div([
    html.Div([
        html.Div([
            html.Label('Gene selection'),
            dcc.Input(
                id='gene-textbox',
                type='text',
                value='',
                placeholder='Select a gene using its ensembl id',
                style={'width': '100%'}
            ),
            dcc.Dropdown(
                id='genes-dropdown',
                options=[ {'label': genes_dropdown_options[key], 'value':key} for key in genes_dropdown_options ],
                placeholder='Select a gene using its name'
            )
        ],
        style={'width': '50%', 'display': 'inline-block'}),

        html.Div([
            html.Img(src='https://frontiersinblog.files.wordpress.com/2018/06/logo.png', style={
                 'height': '100%',
                 'position': 'absolute',
                 'bottom': '-10px',
                 # 'border': '2px solid black',
                 'right': '20px'
            })
        ],
        style={
            'width': '50%',
            'height': '40px',
            'display': 'inline-block',
            # 'border': '2px solid red',
            'position': 'relative',
            'float':'right',
            'bottom': '0'})

    ], style={
        'borderBottom': 'thin lightgrey solid',
        'backgroundColor': 'rgb(250, 250, 250)',
        'padding': '10px 5px'
    }),

    html.Div([

        html.Div([

            dcc.Graph(
                id='tsne-plot',
                figure=cells_scatter_tsne )

        ], style={'width': '50%', 'display': 'inline-block'}),

        html.Div([

            dcc.Graph(
                id='time-strain-plot',
                figure=cells_scatter_plot )

        ], style={'width': '50%', 'display': 'inline-block'})

    ], style={ 'padding': '10px 5px' }),

    html.Div([

        dcc.Graph(
            id='genes-plot',
            figure=genes_scatter_plot )

    ])

],style={'max-width': '1200px', 'margin': '0 auto'})



# trigger on dropdown
@app.callback(
    Output('gene-textbox', 'value'),
    [Input('genes-dropdown', 'value')])
def update_genes_dropdown(dropdown_value):

    if dropdown_value is None:
        return ''
    else:
        return dropdown_value

# trigger on gene plot click
@app.callback(
    Output('genes-dropdown', 'value'),
    [Input('genes-plot', 'clickData')])
def update_genes_dropdown(clickData):

    if clickData is None:
        raise PreventUpdate
    else:
        for option in genes_dropdown_options:
            if genes_dropdown_options[option] == clickData['points'][0]['text']:
                return option

        return clickData['points'][0]['text']

@app.callback(
    Output('time-strain-plot', 'figure'),
    [Input('gene-textbox', 'value')])
def update_time_strain(selected_gene):

    if selected_gene not in genes_options:

        print('default graph')

        cells_scatter_plot = go.Figure(
                data=[],
                layout=go.Layout(
                    # title='US Export of Plastic Scrap',
                    showlegend=True,
                    legend=go.layout.Legend(
                        x=1,
                        y=0
                    ),
                    margin=go.layout.Margin(l=40, r=0, t=40, b=30),
                    xaxis={ 'title': 'strain' },
                    yaxis={ 'title': 'time' }
                )
            )

        for strain in condition['strain'].unique():

            strain_df = condition[ condition['strain'] == strain ]

            # Add first scatter trace with medium sized markers
            cells_scatter_plot.add_trace(
                go.Scatter(
                    mode='markers',
                    x=strain_df['x'],
                    y=strain_df['y'],
                    text= [ v[0] + ', time: ' + str(v[1]) for v in strain_df[['cellname', 'time']].values ],
                    opacity=0.7,
                    marker=dict(
                        # color='Blue',
                        size=10,
                        line={'width': 0.5, 'color': 'white'},
                        symbol=strain_df['shape'].values[0]
                    ),
                    name=strain
                )
            )

        return cells_scatter_plot

    else:

        print('color using selected gene hue')

        gene_reads = np.array(cnt_normalized[ selected_gene ].todense())[0]
        gene_reads = { cell: gene_reads[i] for i, cell in enumerate(cnt_columns) }

        max_read = cnt_normalized[ selected_gene ].max()

        cells_scatter_plot = go.Figure(
                data=[],
                layout=go.Layout(
                    # title='US Export of Plastic Scrap',
                    showlegend=True,
                    legend=go.layout.Legend(
                        x=.87,
                        y=0
                    ),
                    margin=go.layout.Margin(l=40, r=0, t=40, b=30),
                    xaxis={ 'title': 'strain' },
                    yaxis={ 'title': 'time' }
                )
            )

        for strain in condition['strain'].unique():

            strain_df = condition[ condition['strain'] == strain ]

            cell_values = [ gene_reads[cell] for cell in strain_df['cellname'].values ]

            # Add first scatter trace with medium sized markers
            cells_scatter_plot.add_trace(
                go.Scatter(
                    mode='markers',
                    x=strain_df['x'],
                    y=strain_df['y'],
                    text= [ v[0] + ', time: ' + str(v[1]) for v in strain_df[['cellname', 'time']].values ],
                    opacity=0.7,
                    marker=dict(
                        size=10,
                        line={'width': 0.5, 'color': 'white'},
                        symbol=strain_df['shape'].values[0],
                        # set color scale
                        cmax=max_read,
                        cmin=0,
                        color=cell_values,
                        colorbar={ 'title': 'reads' },
                        colorscale="Reds"
                    ),
                    name=strain
                )
            )

        return cells_scatter_plot

@app.callback(
    Output('tsne-plot', 'figure'),
    [Input('gene-textbox', 'value')])
def update_time_strain(selected_gene):

    if selected_gene not in genes_options:
        print('default graph')

        cells_scatter_tsne = go.Figure(
            data=[],
            layout=go.Layout(
                # title='US Export of Plastic Scrap',
                showlegend=True,
                legend=go.layout.Legend(
                    x=1,
                    y=0
                ),
                margin=go.layout.Margin(l=40, r=0, t=40, b=30),
                xaxis={ 'title': 'tSNE 1' },
                yaxis={ 'title': 'tSNE2' }
            )
        )

        for strain in condition['strain'].unique():

            strain_df = condition[ condition['strain'] == strain ]

            # Add first scatter trace with medium sized markers
            cells_scatter_tsne.add_trace(
                go.Scatter(
                    mode='markers',
                    x=strain_df['tSNE_1'],
                    y=strain_df['tSNE_2'],
                    text= [ v[0] + ', time: ' + str(v[1]) for v in strain_df[['cellname', 'time']].values ],
                    opacity=0.7,
                    marker=dict(
                        # color='Blue',
                        size=10,
                        line= {'width': 0.5, 'color': 'white'},
                        symbol=strain_df['shape'].values[0]
                    ),
                    name=strain
                )
            )

        return cells_scatter_tsne

    else:

        print('recolor graph with given hue')

        gene_reads = np.array(cnt_normalized[ selected_gene ].todense())[0]
        gene_reads = { cell: gene_reads[i] for i, cell in enumerate(cnt_columns) }

        max_read = cnt_normalized[ selected_gene ].max()

        cells_scatter_tsne = go.Figure(
            data=[],
            layout=go.Layout(
                # title='US Export of Plastic Scrap',
                showlegend=True,
                legend=go.layout.Legend(
                    x=.87,
                    y=0
                ),
                margin=go.layout.Margin(l=40, r=0, t=40, b=30),
                xaxis={ 'title': 'tSNE 1' },
                yaxis={ 'title': 'tSNE2' }
            )
        )

        for strain in condition['strain'].unique():

            strain_df = condition[ condition['strain'] == strain ]

            cell_values = [ gene_reads[cell] for cell in strain_df['cellname'].values ]

            # Add first scatter trace with medium sized markers
            cells_scatter_tsne.add_trace(
                go.Scatter(
                    mode='markers',
                    x=strain_df['tSNE_1'],
                    y=strain_df['tSNE_2'],
                    text= [ v[0] + ', time: ' + str(v[1]) for v in strain_df[['cellname', 'time']].values ],
                    opacity=0.7,
                    marker=dict(
                        size=10,
                        line= {'width': 0.5, 'color': 'white'},
                        symbol=strain_df['shape'].values[0],
                        # set color scale
                        cmax=max_read,
                        cmin=0,
                        color=cell_values,
                        colorbar={ 'title': 'reads' },
                        colorscale="Reds"
                    ),
                    name=strain
                )
            )

        return cells_scatter_tsne


@app.callback(
    Output('genes-plot', 'figure'),
    [Input('gene-textbox', 'value')])
def update_time_strain(selected_gene):

    if selected_gene not in genes_options:

        print('default graph')

        genes_scatter_plot = go.Figure(
                data=[],
                layout=go.Layout(
                    title='Genes Projection',
                    showlegend=True,
                    legend=go.layout.Legend( x=0, y=1.0 ),
                    margin=go.layout.Margin(l=40, r=0, t=40, b=30),
                    xaxis={ 'title': 'FC strain' },
                    yaxis={ 'title': 'FC time' }
                )
            )

        for organism in strain_time['organism'].unique():

            organism_df = strain_time[ strain_time['organism'] == organism ]

            genes_scatter_plot.add_trace(

                go.Scatter(
                    mode='markers',
                    x=organism_df['fc_strain'],
                    y=organism_df['fc_time'],
                    text=organism_df['Associated.Gene.Name'],
                    opacity=0.8,
                    marker=dict(
                        size=10,
                        line={'width': 0.5, 'color': 'white'},
                        # color = strain_time['organism'],
                        # color='Blue',
                        # symbol=strain_df['shape'].values[0]
                    ),
                    name=organism
                )
            )

        return genes_scatter_plot

    else:

        print('color using selected gene hue')

        import resource

        print('ram consumption', int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000))

        genes_scatter_plot = go.Figure(
                data=[],
                layout=go.Layout(
                    title='Genes Projection',
                    showlegend=True,
                    legend=go.layout.Legend( x=0, y=1.0 ),
                    margin=go.layout.Margin(l=40, r=0, t=40, b=30),
                    xaxis={ 'title': 'FC strain' },
                    yaxis={ 'title': 'FC time' }
                )
            )

        for organism in strain_time['organism'].unique():

            organism_df = strain_time[ strain_time['organism'] == organism ]

            genes_scatter_plot.add_trace(

                go.Scatter(
                    mode='markers',
                    x=organism_df['fc_strain'],
                    y=organism_df['fc_time'],
                    text=organism_df['Associated.Gene.Name'],
                    marker=dict(
                        size=10,
                        line={'width': 0.5, 'color': 'white'},
                        opacity=[ 0.2 if gene_id != selected_gene else 1.0 for gene_id in organism_df['Ensembl.Gene.ID'].values ],
                        # color = strain_time['organism'],
                        # color='Blue',
                        # symbol=strain_df['shape'].values[0]
                    ),
                    name=organism
                )
            )

        return genes_scatter_plot

'''
    dff = df[df['Country Name'] == hoverData['points'][0]['customdata']]
    dff = dff[dff['Indicator Name'] == yaxis_column_name]
    return create_time_series(dff, axis_type, yaxis_column_name)
'''

if __name__ == '__main__':
    app.run_server(debug = True)
    # app.run_server(debug = True, port = 8051)
