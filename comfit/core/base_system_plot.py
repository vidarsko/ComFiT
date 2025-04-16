from typing import Union, Optional, Literal, Tuple, Dict, List, Any, TYPE_CHECKING


import numpy as np
import scipy as sp

from comfit.tool import tool_colormap
from comfit.tool import tool_create_orthonormal_triad
from comfit.tool import tool_set_plot_axis_properties_matplotlib
from comfit.tool import tool_complete_field
from comfit.tool import tool_print_in_color
from comfit.tool import tool_matplotlib_define_3D_plot_ax
from comfit.tool import tool_plotly_define_3D_plot_ax

from comfit.plot import plot_surface_matplotlib
from comfit.plot import plot_field_matplotlib
from comfit.plot import plot_field_plotly
from comfit.plot import plot_complex_field_matplotlib
from comfit.plot import plot_complex_field_plotly
from comfit.plot import plot_angle_field_matplotlib
from comfit.plot import plot_angle_field_plotly
from comfit.plot import plot_vector_field_matplotlib
from comfit.plot import plot_vector_field_plotly
from comfit.plot import plot_field_in_plane_matplotlib
from comfit.plot import plot_field_in_plane_plotly
from comfit.plot import plot_complex_field_in_plane_matplotlib
from comfit.plot import plot_complex_field_in_plane_plotly
from comfit.plot import plot_nodes_matplotlib
from comfit.plot import plot_nodes_plotly
from comfit.plot import plot_vector_field_in_plane_both_plot_libs

from comfit.plot import plot_subplots_matplotlib
from comfit.plot import plot_subplots_plotly

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.tri as mtri
import matplotlib.cm as cm
import matplotlib.axes
import matplotlib.figure

import mpl_toolkits

import plotly.graph_objects as go
import plotly.figure_factory as ff
import plotly.colors as pc

from skimage.measure import marching_cubes

import time

class BaseSystemPlot:
    """Plotting methods for the base system class"""

    def plot_prepare(
        self: 'BaseSystemPlot', 
        field: np.ndarray, 
        field_type: Literal['real', 'complex', 'angle', 'vector'], 
        **kwargs: any
        ) -> Tuple[np.ndarray, Union[go.Figure, plt.Figure], Any, Dict]:
        """Prepare axis and figure for plotting.

        Parameters
        ----------
        field : np.ndarray
            Field to be plotted.
        field_type : str
            Type of the field to be plotted.

        Returns
        -------
        Tuple[np.ndarray, Union[go.Figure, plt.Figure], Any, Dict]
            Tuple containing the field, figure (plotly or matplotlib), 
            axis (plotly or matplotlib) and keyword arguments for plotting.
        """

        plot_lib = kwargs.get('plot_lib', self.plot_lib)
        
        ## Define figure
        if plot_lib == "plotly":
            ax = kwargs.get('ax', {'row': 1, 'col': 1, 'nrows': 1, 'ncols': 1, 'colorbar': False, 'fig': go.Figure()})
            fig = kwargs.get('fig', ax['fig'])

        elif plot_lib == "matplotlib":
            fig = kwargs['fig'] if 'fig' in kwargs else plt.figure()
            ax = kwargs.get('ax', None)

        # Extend field if not a complete array is given
        if field_type in ['real', 'complex', 'angle']:
            field = tool_complete_field(self, field)
        elif field_type == 'vector':
            field_copy = []
            for n in range(len(field)):
                field_copy.append(tool_complete_field(self, field[n]))
            field = np.array(field_copy)

        ## Checks on field
        # Check if the provided field is bool
        if field.dtype == bool:
            field = field.astype(float)

        # Check if the field is nan
        kwargs['field_is_nan'] = False
        if np.all(np.isnan(field)):
            kwargs['field_is_nan'] = True

        # Check if the field is an order parameter containing several components
        if field_type in ['real', 'complex', 'angle']:
            if field.ndim == self.dim+1:
                tool_print_in_color('Warning: The provided field seems to be an order parameter containing several fields. Only the zeroth component will be plotted.')
                field = field[0]

        # Check if the field is complex when it should be real
        if field_type in ['real', 'vector', 'angle']:
            if np.iscomplexobj(field):
                tool_print_in_color("Warning: the provided field was complex. This might be due to residual imaginary parts from the Fourier transform. The imaginary parts will be removed.")
                print('Max imaginary part: ', np.max(np.imag(field)))
                field = np.real(field)

        ## Define default kwarg values
        kwargs['colormap'] = kwargs.get('colormap', 'viridis')
        kwargs['colormap_object'] = tool_colormap(kwargs['colormap'], plot_lib=plot_lib)

        colorbar_default = False
        if field_type == 'complex':
            colorbar_default = True
        if field_type in ['real','angle'] and self.dim>1:
            colorbar_default = True
        if field_type == 'vector' and plot_lib == 'plotly':
            colorbar_default = True
        kwargs['colorbar'] = kwargs.get('colorbar',  colorbar_default)

        kwargs['vlim_symmetric'] = kwargs.get('vlim_symmetric', False)

        grid_default = False
        if field_type == 'real' and self.dim in [1,3]:
            grid_default = True
        if field_type == 'complex':
            grid_default = True
        kwargs['grid'] = kwargs.get('grid', grid_default)

        axis_equal_default = False
        if self.dim in [2,3]:
            axis_equal_default = True
        kwargs['axis_equal'] = kwargs.get('axis_equal', axis_equal_default)

        kwargs['plot_is_3D'] = True if self.dim == 3 else False

        if field_type == 'real':
            # Field limits
            field_min = np.min(field)
            field_max = np.max(field)

            if 'vlim' in kwargs:
                vlim = kwargs['vlim']
                vmin = vlim[0]
                vmax = vlim[1]
            else:
                vmin = field_min
                vmax = field_max
                if 'vlim_symmetric' in kwargs:
                    if kwargs['vlim_symmetric']:
                        vmax = max(abs(vmin), abs(vmax))
                        vmin = -vmax
            
            if plot_lib == 'plotly':
                ax['vmin'] = min(vmin, ax.get('vmin', vmin)) 
                ax['vmax'] = max(vmax, ax.get('vmax', vmax))

            else:
                kwargs['vmin'] = vmin
                kwargs['vmax'] = vmax

        return field, fig, ax, kwargs
    
    def _check_if_fourier_and_adjust(self, field, **kwargs):
        """Check if the field is in Fourier space and adjust it if necessary.

        Parameters
        ----------
        field : np.ndarray
            Field to be checked
            kwargs : Any
        
        Returns
        -------
        Tuple[np.ndarray, Dict]
            Tuple containing the field and the keyword arguments.
        """

        kwargs["fourier"] = kwargs.get("fourier", False)
        if kwargs["fourier"]:
            
            dkx = self.k[0][1]-self.k[0][0]
            phase_shift = 1/(self.xRes*dkx)*np.exp(-1j*self.k[0]*self.xmin)
            if self.dim > 1:
                dky = self.k[1][0,1]-self.k[1][0,0]
                phase_shift = phase_shift*1/(self.yRes*dky)*np.exp(-1j*self.k[1]*self.ymin)
            if self.dim > 2:
                dkz = self.k[2][0,0,1]-self.k[2][0,0,0]
                phase_shift = phase_shift*1/(self.zRes*dkz)*np.exp(-1j*self.k[2]*self.zmin)

            field = np.fft.fftshift(phase_shift*field, axes=range(-self.dim, 0))

            kwargs["x"] = np.fft.fftshift(self.k[0])/(1/self.a0)
            kwargs['xlabel'] = 'kx/(a₀⁻¹)'
            kwargs['xlim'] = kwargs.get('xlim', [np.min(kwargs["x"]), np.max(kwargs["x"])])

            if self.dim > 1:
                kwargs["y"] = np.fft.fftshift(self.k[1])/(1/self.a0)
                kwargs['ylabel'] = 'ky/(a₀⁻¹)'
                kwargs['ylim'] = kwargs.get('ylim', [np.min(kwargs["y"]), np.max(kwargs["y"])])

            if self.dim > 2:
                kwargs["z"] = np.fft.fftshift(self.k[2])/(1/self.a0)
                kwargs['zlabel'] = 'kz/(a₀⁻¹)'
                kwargs['zlim'] = kwargs.get('zlim', [np.min(kwargs["z"]), np.max(kwargs["z"])])

            


        else:

            kwargs['x'] = self.x/self.a0
            kwargs['xlabel'] = kwargs.get('xlabel','x/a₀')
            
            if self.dim > 1:
                kwargs['y'] = self.y/self.a0
                kwargs['ylabel'] = kwargs.get('ylabel','y/a₀')
            if self.dim > 2:
                kwargs['z'] = self.z/self.a0
                kwargs['zlabel'] = kwargs.get('zlabel','z/a₀')

        return field, kwargs

    def plot_field(
            self: 'BaseSystemPlot', 
            field: np.ndarray, 
            **kwargs: Any
            ) -> Tuple[Union[go.Figure, plt.Figure], Any]:
        """Plot a field.

        Parameters
        ----------
        field : np.ndarray
            Field to be plotted.
        kwargs : Any
            Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/ 
            for a full list of keyword arguments.

        Returns
        -------
        Tuple[Union[go.Figure, plt.Figure], Any]
            Tuple containing the figure (plotly or matplotlib) and the axes 
            dictionary (plotly or matplotlib).

        Note
        ----
        plot_lib is a possible keyword argument. If not provided, the default 
        plotting library (plotly or matplotlib) will be used.
        See the corresponding plot_field_plotly or plot_field_matplotlib function
        in the Auxiliary Plot Functions documentation for further details.
        """

        plot_lib = kwargs.get('plot_lib', self.plot_lib)

        field, kwargs = self._check_if_fourier_and_adjust(field, **kwargs)  

        if plot_lib == "plotly":
            return plot_field_plotly(self, field, **kwargs)
        elif plot_lib == "matplotlib":
            return plot_field_matplotlib(self, field, **kwargs)

    def plot_complex_field(
            self: 'BaseSystemPlot', 
            complex_field: np.ndarray, 
            **kwargs: Any
            ) -> Tuple[Union[go.Figure, plt.Figure], Any]:
        """Plot a complex field.

        Parameters
        ----------
        complex_field : np.ndarray
            Complex field to be plotted.
        kwargs : Any
            Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/
            for a full list of keyword arguments.
        
        Returns
        -------
        Tuple[Union[go.Figure, plt.Figure], Any]
            Tuple containing the figure (plotly or matplotlib) and the axes
            dictionary (plotly or matplotlib).
        
        Note
        ----
        plot_lib is a possible keyword argument. If not provided, the default
        plotting library (plotly or matplotlib) will be used.
        See the corresponding plot_complex_field_plotly or plot_complex_field_matplotlib function
        in the Auxiliary Plot Functions documentation for further details.
        """

        plot_lib = kwargs.get('plot_lib', self.plot_lib)

        complex_field, kwargs = self._check_if_fourier_and_adjust(complex_field, **kwargs)  

        if plot_lib == "plotly":
            return plot_complex_field_plotly(self, complex_field, **kwargs)
        elif plot_lib == "matplotlib":
            return plot_complex_field_matplotlib(self, complex_field, **kwargs)


    def plot_angle_field(
            self: 'BaseSystemPlot', 
            angle_field: np.ndarray, 
            **kwargs: Any
            ) -> Tuple[Union[go.Figure, plt.Figure], Any]:
        """Plot an angle field.

        Parameters
        ----------
        angle_field : np.ndarray
            Angle field to be plotted.
        kwargs : Any
            Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/
            for a full list of keyword arguments.
        
        Returns
        -------
        Tuple[Union[go.Figure, plt.Figure], Any]
            Tuple containing the figure (plotly or matplotlib) and the axes
            dictionary (plotly or matplotlib).

        Note
        ----
        plot_lib is a possible keyword argument. If not provided, the default
        plotting library (plotly or matplotlib) will be used.
        See the corresponding plot_angle_field_plotly or plot_angle_field_matplotlib function
        in the Auxiliary Plot Functions documentation for further details.
        """

        plot_lib = kwargs.get('plot_lib', self.plot_lib)

        angle_field, kwargs = self._check_if_fourier_and_adjust(angle_field, **kwargs)
        if kwargs['fourier']:
            raise NotImplementedError("Fourier space angle field plotting not implemented.")

        if plot_lib == "plotly":
            return plot_angle_field_plotly(self, angle_field, **kwargs)
        elif plot_lib == "matplotlib":
            return plot_angle_field_matplotlib(self, angle_field, **kwargs)


    def plot_vector_field(
            self: 'BaseSystemPlot', 
            vector_field: np.ndarray, 
            **kwargs: Any
            ) -> Tuple[Union[go.Figure, plt.Figure], Any]:
        """Plot a vector field.

        Parameters
        ----------
        vector_field : np.ndarray
            Vector field to be plotted.
        kwargs : Any
            Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/
            for a full list of keyword arguments.
        
        Returns
        -------
        Tuple[Union[go.Figure, plt.Figure], Any]
            Tuple containing the figure (plotly or matplotlib) and the axes
            dictionary (plotly or matplotlib).

        Note
        ----
        plot_lib is a possible keyword argument. If not provided, the default
        plotting library (plotly or matplotlib) will be used.
        See the corresponding plot_vector_field_plotly or plot_vector_field_matplotlib function
        in the Auxiliary Plot Functions documentation for further details.
        """

        plot_lib = kwargs.get('plot_lib', self.plot_lib)

        vector_field, kwargs = self._check_if_fourier_and_adjust(vector_field, **kwargs)
        if kwargs['fourier']:
            raise NotImplementedError("Fourier space vector field plotting not implemented.")

        if plot_lib == "plotly":
            return plot_vector_field_plotly(self, vector_field, **kwargs)
        elif plot_lib == "matplotlib":
            return plot_vector_field_matplotlib(self, vector_field, **kwargs)

    def plot_field_in_plane(
        self: 'BaseSystemPlot',
        field: np.ndarray,
        normal_vector: Optional[np.ndarray] = None,
        position: Optional[np.ndarray] = None,
        **kwargs: Any
        ) -> Tuple[Union[go.Figure, plt.Figure], Any]:
        """Plot the field in a plane.

        Parameters
        ----------
        field : np.ndarray
            Field to be plotted.
        normal_vector : np.ndarray, optional
            Normal vector of the plane. If None, the normal vector will be calculated.
        position : np.ndarray, optional
            Position of the plane. If None, the position will be calculated.
        kwargs : Any
            Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/
            for a full list of keyword arguments.
        
        Returns
        -------
        Tuple[Union[go.Figure, plt.Figure], Any]
            Tuple containing the figure (plotly or matplotlib) and the axes
            dictionary (plotly or matplotlib).

        Note
        ----
        plot_lib is a possible keyword argument. If not provided, the default
        plotting library (plotly or matplotlib) will be used.
        See the corresponding plot_field_in_plane_plotly or plot_field_in_plane_matplotlib function
        in the Auxiliary Plot Functions documentation for further details.
        """

        plot_lib = kwargs.get('plot_lib', self.plot_lib)

        field, kwargs = self._check_if_fourier_and_adjust(field, **kwargs)
        if kwargs['fourier']:
            raise NotImplementedError("Fourier space plot in plane not implemented.")

        if plot_lib == "plotly":
            return plot_field_in_plane_plotly(
                self, field, normal_vector, position, **kwargs
            )
        elif plot_lib == "matplotlib":
            return plot_field_in_plane_matplotlib(
                self, field, normal_vector, position, **kwargs
            )

    def plot_complex_field_in_plane(
        self: 'BaseSystemPlot',
        complex_field: np.ndarray,
        normal_vector: Optional[np.ndarray] = None,
        position: Optional[np.ndarray] = None,
        **kwargs: Any
        ) -> Tuple[Union[go.Figure, plt.Figure], Any]:
        """Plot the complex field in a plane.

        Parameters
        ----------
        complex_field : np.ndarray
            Complex field to be plotted.
        normal_vector : np.ndarray, optional
            Normal vector of the plane. If None, the normal vector will be calculated.
        position : np.ndarray, optional
            Position of the plane. If None, the position will be calculated.
        kwargs : Any
            Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/
            for a full list of keyword arguments.

        Returns
        -------
        Tuple[Union[go.Figure, plt.Figure], Any]
            Tuple containing the figure (plotly or matplotlib) and the axes
            dictionary (plotly or matplotlib).

        Note
        ----
        plot_lib is a possible keyword argument. If not provided, the default
        plotting library (plotly or matplotlib) will be used.
        See the corresponding plot_complex_field_in_plane_plotly or plot_complex_field_in_plane_matplotlib function
        in the Auxiliary Plot Functions documentation for further details.
        """

        plot_lib = kwargs.get('plot_lib', self.plot_lib)

        complex_field, kwargs = self._check_if_fourier_and_adjust(complex_field, **kwargs)
        if kwargs['fourier']:
            raise NotImplementedError("Fourier space complex field plot in plane not implemented.")


        if plot_lib == "plotly":
            return plot_complex_field_in_plane_plotly(
                self, complex_field, normal_vector, position, **kwargs
            )
        elif plot_lib == "matplotlib":
            return plot_complex_field_in_plane_matplotlib(
                self, complex_field, normal_vector, position, **kwargs
            )

    def plot_angle_field_in_plane(
        self: 'BaseSystemPlot',
        angle_field: np.ndarray,
        normal_vector: Optional[np.ndarray] = None,
        position: Optional[np.ndarray] = None,
        **kwargs: Any
        ) -> Tuple[Union[go.Figure, plt.Figure], Any]:
        """Plot the angle field in a plane.

        Parameters
        ----------
        angle_field : np.ndarray
            Angle field to be plotted.
        normal_vector : np.ndarray, optional
            Normal vector of the plane. If None, the normal vector will be calculated.
        position : np.ndarray, optional
            Position of the plane. If None, the position will be calculated.
        kwargs : Any
            Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/
            for a full list of keyword arguments.
        
        Returns
        -------
        Tuple[Union[go.Figure, plt.Figure], Any]
            Tuple containing the figure (plotly or matplotlib) and the axes
            dictionary (plotly or matplotlib).

        Note
        ----
        plot_lib is a possible keyword argument. If not provided, the default
        plotting library (plotly or matplotlib) will be used.
        See the corresponding plot_angle_field_plotly or plot_angle_field_matplotlib function
        in the Auxiliary Plot Functions documentation for further details.
        """

        plot_lib = kwargs.get('plot_lib', self.plot_lib)

        angle_field, kwargs = self._check_if_fourier_and_adjust(angle_field, **kwargs)
        if kwargs['fourier']:
            raise NotImplementedError("Fourier space angle field plot in plane not implemented.")


        complex_field = np.exp(1j * angle_field)
        return self.plot_complex_field_in_plane(complex_field, normal_vector=normal_vector, position=position, **kwargs)

    def plot_vector_field_in_plane(
        self: 'BaseSystemPlot',
        vector_field: np.ndarray,
        normal_vector: Optional[np.ndarray] = None,
        position: Optional[np.ndarray] = None,
        spacing = None,
        **kwargs: Any
        ) -> Tuple[Union[go.Figure, plt.Figure], Any]:
        """Plot the vector field in a plane.

        Parameters
        ----------
        vector_field : np.ndarray
            Vector field to be plotted.
        normal_vector : np.ndarray, optional
            Normal vector of the plane. If None, the normal vector will be calculated.
        position : np.ndarray, optional
            Position of the plane. If None, the position will be calculated.
        spacing : float, optional
            Spacing between the vectors.
        kwargs : Any
            Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/
            for a full list of keyword arguments.
        
        Returns
        -------
        Tuple[Union[go.Figure, plt.Figure], Any]
            Tuple containing the figure (plotly or matplotlib) and the axes
            dictionary (plotly or matplotlib).

        Note
        ----
        plot_lib is a possible keyword argument. If not provided, the default
        plotting library (plotly or matplotlib) will be used.
        See the plot_vector_field_both_plot_libs function
        in the Auxiliary Plot Functions documentation for further details.
        """

        plot_lib = kwargs.get('plot_lib', self.plot_lib)

        vector_field, kwargs = self._check_if_fourier_and_adjust(vector_field, **kwargs)
        if kwargs['fourier']:
            raise NotImplementedError("Fourier space vector field plot in plane not implemented.")


        return plot_vector_field_in_plane_both_plot_libs(self, vector_field, normal_vector, position, spacing, **kwargs)

    def plot_nodes(
            self: 'BaseSystemPlot', 
            nodes: Dict[str, Any], 
            **kwargs: Any
            ) ->  Tuple[Union[go.Figure, plt.Figure], Any]:
        """Plot nodes.

        Parameters
        ----------
        nodes : np.ndarray
            Nodes to be plotted.
        kwargs : Any
            Keyword arguments for the plot. See https://comfitlib.com/ClassBaseSystem/
            for a full list of keyword arguments.
        
        Returns
        -------
        Tuple[Union[go.Figure, plt.Figure], Any]
            Tuple containing the figure (plotly or matplotlib) and the axes
            dictionary (plotly or matplotlib).
        
        Note
        ----
        plot_lib is a possible keyword argument. If not provided, the default
        plotting library (plotly or matplotlib) will be used.
        See the corresponding plot_nodes_plotly or plot_nodes_matplotlib function
        in the Auxiliary Plot Functions documentation for further details.
        """

        plot_lib = kwargs.get('plot_lib', self.plot_lib)

        if kwargs.get('fourier', False):
            raise NotImplementedError("Fourier space node does not make sense.")

        if plot_lib == 'plotly':
            return plot_nodes_plotly(self, nodes, **kwargs)
        elif plot_lib == 'matplotlib':
            return plot_nodes_matplotlib(self, nodes, **kwargs)


    # Figure handling methods
    def plot_subplots(
            self: 'BaseSystemPlot', 
            number_of_rows: int, 
            number_of_columns: int, 
            **kwargs: Any
            ) -> Tuple[Union[go.Figure, plt.Figure], Any]:
        """Plot subplots.

        Parameters
        ----------
        number_of_rows : int
            Number of rows in the subplot.
        number_of_columns : int
            Number of columns in the subplot.
        
        Returns
        -------
        Tuple[Union[go.Figure, plt.Figure], Any]
            Tuple containing the figure (plotly or matplotlib) and the axes
            dictionary (plotly or matplotlib).
        """

        plot_lib = kwargs.get('plot_lib', self.plot_lib)

        if plot_lib == "plotly":
            return plot_subplots_plotly(number_of_rows, number_of_columns)
        elif plot_lib == "matplotlib":
            return plot_subplots_matplotlib(number_of_rows, number_of_columns, **kwargs)

    def plot_save(
        self: 'BaseSystemPlot', 
        fig: Union[plt.Figure, go.Figure], 
        counter: Optional[int] = None, 
        **kwargs: Any
        ) -> None:
        """Save a figure.

        Parameters
        ----------
        fig : Union[plt.Figure, go.Figure]
            Figure to be saved.
        counter : int, optional
            Counter for the figure. If None, the figure will be saved as 'plot.png'.
        kwargs : Any
            Optional arguments: ID, (matplotlib:) image_size_inches, dpi. 
            (plotly:) width, height.
        
        Returns
        -------
        None
            Save the figure. 
        """

        # Check if the figure is the first argument (deprecated)
        if isinstance(fig, int):
            tool_print_in_color("Warning: The plot_save method has been called with the figure as the second argument. This is deprecated. Please call the method with the figure as the first argument, then counter.")
            counter_tmp = fig
            fig = counter
            counter = counter_tmp
        
        # Get ID
        ID=kwargs.get('ID', None)

        # Keyword arguments
        image_size_inches=kwargs.get('image_size_inches', (6,5))
        dpi=kwargs.get('dpi', 100)

        # Set filename
        if counter is None:
            if ID is None:
                filename = 'plot.png'
            else:
                filename = f'plot_{ID}.png'
        else:
            if ID is None:
                filename = f'plot_{counter}.png'
            else:
                filename = f'plot_{counter}_{ID}.png'


        # Save the figure
        # Determine the plotting library based on figure type
        if isinstance(fig, plt.Figure):
            plot_lib = "matplotlib"
        elif isinstance(fig, go.Figure):
            plot_lib = "plotly"
        else:
            plot_lib = self.plot_lib

        if plot_lib == "plotly":
            width=kwargs.get('width', 800)
            height=kwargs.get('height', 600)
            fig.write_image(filename, width=width, height=height)

        elif plot_lib == "matplotlib":
            fig.set_size_inches(image_size_inches)
            fig.savefig(filename, dpi=dpi)
            plt.close(fig)


    def show(
            self: 'BaseSystemPlot', 
            fig: Union[plt.Figure, go.Figure]
            ) -> None:
        """Show a figure.

        Parameters
        ----------
        fig : Union[plt.Figure, go.Figure]
            Figure to be shown.
        
        Returns
        -------
        None
            Show the figure.
        """

        plot_lib = "plotly" if isinstance(fig, go.Figure) else "matplotlib"

        if plot_lib == "matplotlib":
            plt.show()
        elif plot_lib == "plotly":
            fig.show()