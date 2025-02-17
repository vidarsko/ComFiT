import numpy as np
import re

import plotly.graph_objects as go

def tool_plotly_find_next_xN(fig):
    """
    Find the next xaxis in a set of axs
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