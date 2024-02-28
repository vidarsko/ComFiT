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
| `vlim` | List or tuple consisting of the lower and upper limit of the value to be plotted. Only relevant for `plot_field`. | None |
| `vlim_symmetric` | A Boolean parameter specifying whether the value limits should be symmetric. Only relevant for `plot_field`. | `False` |
| `colorbar` | Boolean parameter indicating whether or not to plot the colorbar | `True` (if applicable)|
| `colormap` | String specifying the colormap to be used | Varies |
| `grid` | Boolean parameter indicating whether or not to plot the axes grid | `False` |
| `hold` | Boolean parameter indicating whether or not to hold the current plot | `False` |
| `plot_shadows` | Boolean parameter indicating whether or not to plot the shadows of the objects. Only applicable for `plot_complex_field`. | `True` |
| `fig` | `matplotlib` figure handle | None|
| `ax` | `matplotlib` axis handle | None|
| `xticks` | List of ticks on the x-axis | None |
| `xticklabels` | List of labels for the ticks on the x-axis | None |
| `yticks` | List of ticks on the y-axis | None |
| `yticklabels` | List of labels for the ticks on the y-axis | None |
| `zticks` | List of ticks on the z-axis | None |
| `zticklabels` | List of labels for the ticks on the z-axis | None |


## Plotting functions

The `BaseSystem` class comes pre-packaged with a number of plotting functions to plot four different types of fields.

### Real fields

Real fields are fields that take real values, for example the temperature in a room.

#### `plot_field`

The `plot_field` function is used to plot a real field.

??? note "Example"

    ```python
    import comfit as cf
    import matplotlib.pyplot as plt
    import numpy as np

    fig = plt.figure()

    ax1 = fig.add_subplot(131)
    bs = cf.BaseSystem(1,xRes=31)
    field = bs.x**2
    bs.plot_field(field,ax=ax1)

    ax2 = fig.add_subplot(132)
    bs = cf.BaseSystem(2,xRes=31,yRes=31)
    field = bs.x**2 + bs.y**2
    bs.plot_field(field,ax=ax2)

    ax3 = fig.add_subplot(133, projection='3d')
    bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
    field = bs.x**2 + bs.y**2 + bs.z**2
    bs.plot_field(field,ax=ax3)

    plt.show()
    ```

    ![](images/plotting_plot_field_demo.png#only-light)
    ![](images/plotting_plot_field_demo-colorinverted.png#only-dark)

#### `plot_field_in_plane` 

The `plot_field_in_plane` function is used to plot a real field in a plane.

??? note "Example"
    ```python 
    import comfit as cf
    import matplotlib.pyplot as plt
    import numpy as np


    fig = plt.figure()

    ax1 = fig.add_subplot(121, projection='3d')
    bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
    field = (bs.x**2 + bs.y**2 + bs.z**2)
    bs.plot_field_in_plane(field, ax=ax1)

    ax2 = fig.add_subplot(122, projection='3d')
    bs.plot_field_in_plane(field, ax=ax2, normal_vector=[1,1,0],position=[10,10,10])

    plt.show()
    ```

    ![](images/plotting_plot_field_in_plane_demo.png#only-light)
    ![](images/plotting_plot_field_in_plane_demo-colorinverted.png#only-dark)

### Complex fields

Complex fields are fields that take complex values, for example the electric field in a light wave.

#### `plot_complex_field` 

The `plot_complex_field` function is used to plot a complex field.

??? note "Example"

    ```python
    import comfit as cf
    import matplotlib.pyplot as plt
    import numpy as np

    fig = plt.figure()

    ax1 = fig.add_subplot(231)
    bs = cf.BaseSystem(1,xRes=31)
    field = bs.x**2*np.exp(1j*bs.x/3)
    bs.plot_complex_field(field,ax=ax1)

    ax2 = fig.add_subplot(232)
    bs = cf.BaseSystem(2,xRes=31,yRes=31)
    field = (bs.x**2 + bs.y**2)*np.exp(1j*bs.x/3)
    bs.plot_complex_field(field,ax=ax2,plot_method='phase_angle')

    ax3 = fig.add_subplot(233, projection='3d')
    bs = cf.BaseSystem(2,xRes=31,yRes=31)
    field = (bs.x**2 + bs.y**2)*np.exp(1j*bs.x/3)
    bs.plot_complex_field(field,ax=ax3,plot_method='3Dsurface')

    ax5 = fig.add_subplot(235, projection='3d')
    bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
    field = (bs.x**2 + bs.y**2 + bs.z**2)*np.exp(1j*bs.x/3)
    bs.plot_complex_field(field,ax=ax5,plot_method='phase_angle')

    ax6 = fig.add_subplot(236, projection='3d')
    bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
    field = (bs.x**2 + bs.y**2 + bs.z**2)*np.exp(1j*bs.x/3)
    bs.plot_complex_field(field,ax=ax6,plot_method='phase_blob')

    plt.show()
    ```

    ![](images/plotting_plot_complex_field_demo.png#only-light)
    ![](images/plotting_plot_complex_field_demo-colorinverted.png#only-dark)

#### `plot_complex_field_in_plane` 

The `plot_complex_field_in_plane` function is used to plot a complex field in a plane.
The modulus of the complex field is shown as the alpha channel, where the minimum modulus value is transparent and the maximum modulus value is opaque.
The phase of the complex field is shown as the color of the field, where the color is determined by the angle color scheme.

??? note "Example"
    ```python 
    import comfit as cf
    import matplotlib.pyplot as plt
    import numpy as np

    fig = plt.figure()

    ax1 = fig.add_subplot(121, projection='3d')
    bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
    complex_field = (bs.x**2 + bs.y**2 + bs.z**2)*np.exp(1j*bs.y/3)
    bs.plot_complex_field_in_plane(complex_field, ax=ax1)

    ax2 = fig.add_subplot(122, projection='3d')
    bs.plot_complex_field_in_plane(complex_field, ax=ax2, normal_vector=[0,0,1],position=[10,10,10])

    plt.show()
    ```

    ![](images/plotting_plot_complex_field_in_plane_demo.png#only-light)
    ![](images/plotting_plot_complex_field_in_plane_demo-colorinverted.png#only-dark)

### Angle fields

Angle fields are fields that take values in the interval $[-\pi,\pi]$, for example the phase of a complex field.

#### `plot_angle_field`

The `plot_angle_field` function is used to plot an angle field.

??? note "Example"
    ```python
    import comfit as cf
    import matplotlib.pyplot as plt
    import numpy as np

    fig = plt.figure()

    ax1 = fig.add_subplot(131)
    bs = cf.BaseSystem(1,xRes=31)
    angle_field = np.mod((bs.x)/5,2*np.pi)-np.pi
    bs.plot_angle_field(angle_field,ax=ax1)

    ax2 = fig.add_subplot(132)
    bs = cf.BaseSystem(2,xRes=31,yRes=31)
    angle_field = np.mod((bs.x + 2*bs.y)/5,2*np.pi)-np.pi
    bs.plot_angle_field(angle_field,ax=ax2)

    ax3 = fig.add_subplot(133, projection='3d')
    bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
    angle_field = np.mod((bs.x + 2*bs.y + 3*bs.z)/5,2*np.pi)-np.pi
    bs.plot_angle_field(angle_field,ax=ax3)

    plt.show()
    ```

    ![](images/plotting_plot_angle_field_demo.png#only-light)
    ![](images/plotting_plot_angle_field_demo-colorinverted.png#only-dark)

#### `plot_angle_field_in_plane` 

The `plot_angle_field_in_plane` function is used to plot an angle field in a plane.

??? note "Example"
    ```python 
    import comfit as cf
    import matplotlib.pyplot as plt
    import numpy as np

    fig = plt.figure()

    ax1 = fig.add_subplot(121, projection='3d')
    bs = cf.BaseSystem(3,xRes=31,yRes=31,zRes=31)
    angle_field = np.mod((bs.x + 2*bs.y + 3*bs.z)/5,2*np.pi)-np.pi
    bs.plot_angle_field_in_plane(angle_field, ax=ax1)

    ax2 = fig.add_subplot(122, projection='3d')
    bs.plot_angle_field_in_plane(angle_field, ax=ax2, normal_vector=[0,0,1],position=[10,10,10])

    plt.show()
    ```

    ![](images/plotting_plot_angle_field_in_plane_demo.png#only-light)
    ![](images/plotting_plot_angle_field_in_plane_demo-colorinverted.png#only-dark)

### Vector fields

#### `plot_vector_field`

The `plot_vector_field` function is used to plot a vector field.
Vector fields are usually plotted blue.
Together with the typical keyword arguments, the `plot_vector_field` function has the kewyword `spacing` which determines the spacing between the arrows in the plot.

The behavior of this plot function is dependent on the interplay between the dimension of the system and the dimension $n$ of the vector field.

| System dimension | $n=1$ | $n=2$ | $n=3$ |
| ----------------- | ----- | ----- | ----- |
|`dim=1` |  Plotted as a 1D plot. | Plotted along the x-axis with the x-component of the vector field as the y-component of the plot and the y-component of the vector field as the z-component of the plot. |  Plotted along the x-axis. The x- and y- and z- component of the vector field are plotted along the x-, y- and z-component. The y- and z- component are shown as is in the perpendicular dimensions, while the x-component has been normalized by the max norm of the input vector field and subsequently scaled by `xmax/3`. The funciton is made for visualization purposes in mind, and quantitative analysis should plot each component individually. |


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