from plotly.subplots import make_subplots

def tool_make_subplots(number_of_rows, number_of_columns, *args):
    """Makes subplots from a series of plotly figures.

    Args:
        *args: A series of plotly figures.
    Returns:
        A plotly figure with subplots.
    """
    if len(args) != number_of_rows * number_of_columns:
        raise Exception("The number of subplots does not match the number of rows and columns.")

    specs_list = [[{} for _ in range(number_of_columns)] for _ in range(number_of_rows)]

    #List of known 3D plot types in Plotly
    types_of_3D_traces = {'scatter3d', 'surface', 'mesh3d', 'cone', 'streamtube', 'isosurface', 'volume'}

    n = 0
    for row in specs_list:
        for spec in row:
            subplot_fig = args[n]
            if len(subplot_fig.data) > 1:
                print('\033[91m Warning: Multiple traces in the same subplot. This might cause issues. \033[0m')

            for trace in subplot_fig.data: #MIght be some issue with different types of traces in the same figure (Vidar 12.09.24)
                if trace.type in types_of_3D_traces:
                    spec['type'] = 'surface'
                else:
                    spec['type'] = 'xy'
            n += 1

    fig = make_subplots(rows=number_of_rows, cols=number_of_columns, specs=specs_list)

    n = 0
    for row in range(1, number_of_rows + 1):
        for col in range(1, number_of_columns + 1):
            subplot_fig = args[n]
            for trace in subplot_fig.data:
                fig.add_trace(trace, row=row, col=col)
                
            # fig.update_layout(subplot_fig.layout)

            # Update x and y axes layout for each subplot (for 2D plots)
            fig.update_xaxes(subplot_fig.layout.xaxis, row=row, col=col)
            fig.update_yaxes(subplot_fig.layout.yaxis, row=row, col=col)

            # Set aspect ratio for 2D plots (if required)
            if all(trace.type not in types_of_3D_traces for trace in subplot_fig.data):
                if hasattr(subplot_fig.layout, 'yaxis') and 'scaleanchor' in subplot_fig.layout.yaxis:
                    if subplot_fig.layout.yaxis.scaleanchor == 'x':
                        fig.update_yaxes(scaleanchor='x', row=row, col=col)
            
            # Update 3D scene layout for each subplot (for 3D plots)
            if any(trace.type in types_of_3D_traces for trace in subplot_fig.data):
                fig.update_scenes(subplot_fig.layout.scene, row=row, col=col)

            n += 1

    

    return fig

