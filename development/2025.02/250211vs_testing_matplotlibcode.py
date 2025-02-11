import matplotlib.pyplot as plt
import numpy as np

def create_subplot_grid(figures, rows, cols):
    """
    Arrange multiple matplotlib figures into a subplot grid.
    
    Args:
        figures (list): List of matplotlib Figure objects
        rows (int): Number of rows in subplot grid
        cols (int): Number of columns in subplot grid
    
    Returns:
        matplotlib.figure.Figure: Combined figure with subplots
    """
    combined_fig = plt.figure(figsize=(cols * 6, rows * 4))
    
    for idx, fig in enumerate(figures):
        if idx >= rows * cols:
            break
            
        # Create subplot
        ax = combined_fig.add_subplot(rows, cols, idx + 1)
        
        # Get the axes content from original figure
        original_ax = fig.get_axes()[0]
        
        # Copy line plots
        for line in original_ax.get_lines():
            ax.plot(line.get_xdata(), line.get_ydata(), 
                   color=line.get_color(),
                   linestyle=line.get_linestyle(),
                   marker=line.get_marker())
            
        # Copy scatter plots
        for collection in original_ax.collections:
            if isinstance(collection, plt.matplotlib.collections.PathCollection):
                ax.scatter(collection.get_offsets()[:, 0], 
                          collection.get_offsets()[:, 1],
                          c=collection.get_facecolor())
        
        # Copy axis labels and title
        ax.set_xlabel(original_ax.get_xlabel())
        ax.set_ylabel(original_ax.get_ylabel())
        ax.set_title(original_ax.get_title())
        
        # Copy axis limits
        ax.set_xlim(original_ax.get_xlim())
        ax.set_ylim(original_ax.get_ylim())
    
    combined_fig.tight_layout()
    return combined_fig



fig1 = plt.figure()
plt.plot([1, 2, 3], [4, 5, 6])

fig2 = plt.figure()
plt.scatter([1, 2, 3], [6, 5, 4])

combined = create_subplot_grid([fig1, fig2], rows=1, cols=2)
plt.show()