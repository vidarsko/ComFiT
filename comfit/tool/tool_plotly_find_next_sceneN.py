# General packages
import plotly.graph_objects as go


def tool_plotly_find_next_sceneN(
        fig : go.Figure
        ) -> str:
    """Find the next available scene name in a Plotly figure.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to analyze.

    Returns
    -------
    str
        The next available scene name. Returns 'scene' if no scenes are used,
        otherwise returns 'sceneN' where N is one more than the current highest scene number.
    """
    
    used_scene = set()
    
    # Check all traces in the figure
    for trace in fig.data:
        if hasattr(trace, 'scene'):
            used_scene.add(trace.scene)
    
    return 'scene' if not used_scene else f'scene{len(used_scene)+1}'


if __name__ == "__main__":
    import plotly.graph_objects as go
    
    fig = go.Figure()
    fig.add_trace(go.Scatter3d(x=[1, 2, 3], y=[4, 5, 6], z=[7, 8, 9], scene='scene'))
    fig.add_trace(go.Scatter3d(x=[2, 3, 4], y=[5, 6, 7], z=[8, 9, 10], scene='scene2'))

    print(tool_plotly_find_next_sceneN(fig))
    print('Expected: scene3')

    fig = go.Figure()
    print(tool_plotly_find_next_sceneN(fig))
    print('Expected: scene')