
def tool_plotly_find_next_sceneN(fig):
    """
    Find the next scene in a set of scenes
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