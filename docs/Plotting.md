# Plotting

The package comes with a lot of plotting functions, so a general exposition is a good idea.

The default plot command is simply `plot()`, which will plot the current
state of the system according to some arbitrarily chosen standard.

While the former is more widely used and documented, the latter is better for 3D visualizations.
The libraries have different nomenclature for the objects that go into the plot functions, which is useful to learn in order to make the behaviour as expected.

Both `matplotlib` uses `figure` to designate the interface on which the plot is drawn.
A new figure is produced by either of the following commands

```python
import matplotlib.pyplot as plt

fig1 = plt.fig()
```
In matplotlib, one level under, we find `axes` handle.
A `matplotlib` figure can contain multiply different `axes` as in several subplots.
The `axes` object is made as follows

```python
ax = fig.add_subplot(111)
```

If you are plotting a 3D object, then you will need to specify that

```python
ax = fig.add_subplot(111, projection='3d')
```

which will construct a 3D `axes` object.

The standard is that when a plotting function is called *without* a keyword argument specifying the current figure or axes, then the current figure will be cleared and potential axes (in the case of matplotlib) will be created onto it.

If a figure is provided by the keyword `fig` with matplotlib, then it will be cleared and the new plot will be plotted on it.
This is because with no reference to which axes the plot is meant to go ontop, there is no way of knowing.

If an axes object is provided by the keyword `ax`, then the new plot will be plotted on top of that axes object if not the keyword argument `hold=False` is also provided.

To show the current plot, one writes

```python
plt.show()
```

which will pause the simulation untill the plot window has been closed.
In order to draw the image and continue the simulation, one needs to write

```python
plt.draw()
plt.pause(0.01)
```


## Plotting keywords

The following list gives the keyword arguments that determine the layout of the resulting plot.
These keywords can be passed to any plot function.
`bs` refers to an instance of the `BaseSystem` class.

| Keyword         | Definition         | Default value |
| ------------------ | --------------- | ----------- |
| `xlabel` | The label on the x-axis | $ x/a_0$|
| `ylabel` | The label on the y-axis | $d=1$: None |
| | | $d = 2$: $y/a_0$|
| | | $d = 3$: $y/a_0$|
| `zlabel` | The label on the z-axis | $d = 1$: None |
| | | $d=2$: None |
| | | $d=3$: $z/a_0$ |
| `suptitle` | The figure title | None |
| `title` | The axes title | None|
| `xmin` | The lower limit on the x-axis | `bs.xmin` |
| `xmax`| The upper limit on the x-axis | `bs.xmax - bs.dx` |
| `xlim`| A list or tuple consisting of the lower and upper limit on the x-axis. If `xlim` is provided, it trumps any provided `xmin` or `xmax`. | None|
| `ymin` | The lower limit on the x-axis | $d=1$:  None |
| | | $d = 2$: `bs.ymin` |
| | | $d = 3$: `bs.ymin` |
| `ymax`| The upper limit on the x-axis | $d=1$: None |
| | | $d = 2$: `bs.ymax-bs.dy` |
| | | $d = 2$: `bs.ymax-bs.dy` |
| `ylim`| A list or tuple consisting of the lower and upper limit on the y-axis. If `ylim` is provided, it trumps any provided `ymin` or `ymax`. | None |
| `zmin` | The lower limit on the x-axis | $d=1$:  None |
| | | $d = 2$: None |
| | | $d = 3$: `bs.zmin` |
| `zmax`| The upper limit on the x-axis | $d=1$: None |
| | | $d = 2$: None |
| | | $d = 2$: `bs.zmax-bs.dz` |
| `zlim`| List or tuple consisting of the lower and upper limit on the z-axis. If `zlim` is provided, it trumps any provided `zmin` or `zmax`. | None |
| `vmin` | Lower limit on the field to be plotted. In the case of a complex function, this is the lower limit of the absolute value of the field to be plotted. |None|
| `vmax` | Upper limit on the value of field to be plotted. In the case of a complex function, this is the upper limit of the absolute value of the field to be plotted. |None|
| `vlim` | List or tuple consisting of the lower and upper limit of the value to be plotted. | None |
| `vlim_symmetric` | A Boolean parameter specifying whether the value limits should be symmetric | `False` |
| `colorbar` | Boolean parameter indicating whether or not to plot the colorbar | `True` (if applicable)|
| `colormap` | String specifying the colormap to be used | Varies |
| `grid` | Boolean parameter indicating whether or not to plot the axes grid | `False` |
| `plotting_lib` | String specifying the plotting library to be used for visualization. | `matplotlib` |
| `ax` | `matplotlib` axis handle | None|

## Predefined figure sizes

Obviously, a lot of these plots are meant for publication. Therefore,
there is a tool package that contains predefined and fitting sizes for
the figures in question. $$TO BE ADDED$$

## Animation

## Angle color scheme

In many of the plotting functions, we are plotting angles, for example in plotting the phase
of a complex number or the value of an order parameter on S1
. In these cases, all values
modulus 2π are eqvuivalent, but if one uses a regular color scheme, this equivalence is not
readily visible. Therefore, when expressing angles, we use the color scheme shown in Fig. 1.1.
This has the benefit of wrapping around itself at θ = ±π, stressing that these correspond

![Angle color scheme](images/conventions_angle_colormap.png#only-light)
![Angle color scheme](images/conventions_angle_colormap-colorinverted.png#only-dark)

*Angle color scheme.* The color scheme follows the hsv color circle going through  $\theta=0$ (Red), $\theta=\pi/3$ (Yellow), $\theta=2\pi/3$ (Lime), $\theta = \pm \pi$ (Aqua), $\theta = -2\pi/3$ (Blue), $\theta = -\pi/3$ (Fuchsia).