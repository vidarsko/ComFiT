from comfit.tool import tool_plotly_find_next_xN

def tool_plotly_define_2D_plot_ax(ax, fig): #Defines xN, yN and plot_dimension
    #Plot nature
    ax['xN'] = ax.get('xN', tool_plotly_find_next_xN(fig))
    ax['yN'] = 'y' if ax['xN'] == 'x' else f'y{ax["xN"][1:]}'

    # Plot dimension
    ax['plot_dimension'] = ax.get('plot_dimension', 2)
    if ax['plot_dimension'] != 2:
        raise ValueError(f"Plot dimension must be 2 for 2D plots. Got {ax['plot_dimension']}")

    # Axis number
    axis_number = 1 if ax['xN'] == 'x' else int(ax['xN'][1:])  # extracts everything after 'x' and converts to integer
    
    ax['xaxisN'] = f'xaxis{axis_number}'
    ax['yaxisN'] = f'yaxis{axis_number}'

    return ax