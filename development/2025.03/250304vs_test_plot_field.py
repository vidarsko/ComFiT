import comfit as cf

# Initialize BaseSystem object
bs = cf.BaseSystem(2, plot_lib='plotly')

colormap = 'sunburst'

# Create field
field = (bs.x+bs.y)/200

bs.plot_lib = 'plotly'
fig, ax = bs.plot_field(field, colorbar=True, colormap = colormap)

bs.plot_lib = 'matplotlib'
fig, ax = bs.plot_field(field, colorbar=True, colormap = colormap)

bs.show(fig)