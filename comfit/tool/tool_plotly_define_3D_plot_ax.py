from comfit.tool.tool_plotly_find_next_sceneN import tool_plotly_find_next_sceneN

def tool_plotly_define_3D_plot_ax(ax,fig):

    ax['plot_dimension'] = ax.get('plot_dimension', 3)
    if ax['plot_dimension'] != 3:
        raise ValueError(f"Plot dimension must be 3 for 3D plots. Got {ax['plot_dimension']}")

    ax['sceneN'] = ax.get('sceneN', tool_plotly_find_next_sceneN(fig))

    return ax