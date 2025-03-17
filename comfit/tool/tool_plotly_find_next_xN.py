# General packages 
import numpy as np
import plotly.graph_objects as go

def tool_plotly_find_next_xN(
        fig : go.Figure
        ) -> str:
    """Find the next available x-axis identifier in a Plotly figure.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to analyze.

    Returns
    -------
    str
        The next available x-axis identifier. Returns 'x' if no x-axes are used,
        or 'xN' where N is one more than the current highest axis number.
    """
    used_xaxis = set()
    
    # Check all traces in the figure
    for trace in fig.data:
        if hasattr(trace, 'xaxis'):
            used_xaxis.add(trace.xaxis)
    
    return 'x' if not used_xaxis else f'x{len(used_xaxis)+1}'

if __name__ == "__main__":
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[1, 2, 3], y=[4, 5, 6], xaxis='x1'))
    fig.add_trace(go.Scatter(x=[1, 2, 3], y=[4, 5, 6], xaxis='x2'))
    fig.add_trace(go.Scatter(x=[1, 2, 3], y=[4, 5, 6], xaxis='x3'))
    fig.add_trace(go.Scatter(x=[1, 2, 3], y=[4, 5, 6], xaxis='x4'))
    fig.add_trace(go.Scatter(x=[1, 2, 3], y=[4, 5, 6], xaxis='x5'))
    fig.add_trace(go.Scatter(x=[1, 2, 3], y=[4, 5, 6], xaxis='x6'))
    print(tool_plotly_find_next_xN(fig)) 
    print('Expected: x7')

    fig = go.Figure()
    print(tool_plotly_find_next_xN(fig)) 
    print('Expected: x')