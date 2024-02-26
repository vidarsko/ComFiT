# Plotting

The `ComFiT` package uses `matplotlib` as the default plotting library.
This is because `matplotlib` is a very versatile and well-documented library that is widely used in the scientific community.
In order to master the use of `matplotlib`, one needs to understand the basic structure of the library.

The basic structure of `matplotlib` is that it has a `figure` object that contains `axes` objects.
One can think of the `figure` object as the window in which the plot is drawn, and the `axes` object as the plot itself.
A new figure is made as follows

```python
import matplotlib.pyplot as plt

fig1 = plt.figure()
```

Given a figure object, one may define a new `axes` object on it as follows

```python
ax = fig.add_subplot(111)
```
The `111` is a shorthand for `1,1,1` and means that the 
If you are plotting a 3D object, then you will need to specify that

```python
ax = fig.add_subplot(111, projection='3d')
```

which will construct a 3D `axes` object.

The convention followed in `ComFiT` are as follows:

* When a plotting function is called *without* a keyword argument specifying the current figure or axes, then the current figure will be cleared and potential axes (in the case of matplotlib) will be created onto it.
This is because with no reference to which axes the plot is meant to go ontop, there is no way of knowing.
* If a figure is provided by the keyword `fig=myfig` with, then it will be cleared and the new plot will be plotted on `myfig`.
This is because with no reference to which axes the plot is meant to go ontop, there is no way of knowing.
* If an axes object is provided by the keyword `ax`, then the `ax` instance will be cleared and the new plot will be plotted on `ax`, unless the keyword `hold=True` is provided, in which case the new plot will be plotted ontop of the old plot.

To show the current plot, one writes

```python
plt.show()
```

which will pause the simulation until the plot window has been closed.
In order to draw the image and continue the simulation, as for instance when viewing a simulation live, one needs to write

```python
plt.pause(0.01)
```


## Plotting keywords

The following list gives the keyword arguments that determine the layout of the resulting plot.
These keywords can be passed to any plot function.
`bs` refers to an instance of the `BaseSystem` class.
In some cases, default values of other parameter depend on the value of `dim`, and are represented by curly brackets:

$$
\left \lbrace \begin{array}{l} \textrm{default value if } \texttt{dim }= 1 \\ \textrm{default value if } \texttt{dim }= 2  \\ \textrm{default value if } \texttt{dim }= 3  \\ \end{array} \right \rbrace
$$

| Keyword         | Definition         | Default value |
| ------------------ | --------------- | ----------- |
| `xlabel` | The label on the x-axis | $x/a_0$|
| `ylabel` | The label on the y-axis |  $\left \lbrace \begin{array}{c} \texttt{none} \\ y/a_0 \\  y/a_0 \\ \end{array} \right \rbrace$  |
| `zlabel` | The label on the z-axis |  $\left \lbrace \begin{array}{c} \texttt{none} \\ \texttt{none} \\  z/a_0 \\ \end{array} \right \rbrace$  |
| `suptitle` | The figure title | None |
| `title` | The axes title | None|
| `xmin` | The lower limit on the x-axis | `bs.xmin` |
| `xmax`| The upper limit on the x-axis | `bs.xmax - bs.dx` |
| `xlim`| A list or tuple consisting of the lower and upper limit on the x-axis. If `xlim` is provided, it trumps any provided `xmin` or `xmax`. | None|
| `ymin` | The lower limit on the y-axis | $\left \lbrace \begin{array}{c} \texttt{none} \\ \texttt{bs.ymin} \\  \texttt{bs.ymin} \\ \end{array} \right \rbrace$ |
| `ymax`| The upper limit on the y-axis | $\left \lbrace \begin{array}{c} \texttt{none} \\ \texttt{bs.ymax-bs.dy} \\  \texttt{bs.ymax-bs.dy} \\ \end{array} \right \rbrace$ |
| `ylim`| A list or tuple consisting of the lower and upper limit on the y-axis. If `ylim` is provided, it trumps any provided `ymin` or `ymax`. | None |
| `zmin` | The lower limit on the z-axis | $\left \lbrace \begin{array}{c} \texttt{none} \\ \texttt{none} \\  \texttt{bs.zmin} \\ \end{array} \right \rbrace$ |
| `zmax`| The upper limit on the z-axis | $\left \lbrace \begin{array}{c} \texttt{none} \\ \texttt{none} \\  \texttt{bs.zmax-bs.dz} \\ \end{array} \right \rbrace$ |
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