import comfit as cf


bs = cf.BaseSystem(2, plot_lib='plotly', xlim=[-10,10], ylim=[-10,10])
colormap = 'sunburst'
# Create field
field = (bs.x+bs.y)
fig, ax = bs.plot_angle_field(field, colorbar=True, colormap=colormap, title=r'$\theta$')
bs.show(fig)
# self.assertIsNotNone(fig, msg=f'plot_angle_field in 2D with custom colorbar failed for plot_lib={plot_lib}')