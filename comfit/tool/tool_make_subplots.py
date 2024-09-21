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
    n_2d = -1 # Counter for 2D plots
    n_3d = -1 # Counter for 3D plots
    for row in range(1, number_of_rows + 1):
        for col in range(1, number_of_columns + 1):
            subplot_fig = args[n]
            plot_is_3D = any(trace.type in types_of_3D_traces for trace in subplot_fig.data)

            if plot_is_3D:
                n_3d += 1
            else:
                n_2d += 1
                xaxis_id = f'x{n_2d+1}'
                yaxis_id = f'y{n_2d+1}'


            for trace in subplot_fig.data:
                if hasattr(trace, 'colorbar') and trace.colorbar is not None:
                    x_offset=0.15
                    trace.colorbar.x = x_offset + (1-x_offset)*col/number_of_columns  # Adjust position based on column
                    y_offset=0.3
                    trace.colorbar.y = y_offset + (1-y_offset)*(1-row/number_of_rows)  # Centered vertically
                    # trace.colorbar.yanchor = 'middle'  # Ensure the colorbar is centered
                fig.add_trace(trace, row=row, col=col)
                

            # Update x and y axes layout for each subplot (for 2D plots)
            if not plot_is_3D:
                fig.update_xaxes(subplot_fig.layout.xaxis, row=row, col=col)
                fig.update_yaxes(subplot_fig.layout.yaxis, row=row, col=col)
                fig.update_yaxes(scaleanchor=xaxis_id, scaleratio=1, row=row, col=col)
                # fig.update_xaxes(scaleanchor=yaxis_id, scaleratio=1, row=row, col=col)
                # Explicitly enable gridlines
                fig.update_xaxes(showgrid=True, row=row, col=col)
                fig.update_yaxes(showgrid=True, row=row, col=col)  

                
            
            # Update 3D scene layout for each subplot (for 3D plots)
            if plot_is_3D:
                fig.update_scenes(subplot_fig.layout.scene, row=row, col=col)

            n += 1

    return fig

