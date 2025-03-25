from comfit.core.base_system_init import BaseSystemInit
from comfit.core.base_system_conf import BaseSystemConf
from comfit.core.base_system_evolve import BaseSystemEvolve
from comfit.core.base_system_calc import BaseSystemCalc
from comfit.core.base_system_plot import BaseSystemPlot
from comfit.core.base_system_get import BaseSystemGet

import scipy as sp

class BaseSystem(BaseSystemInit, BaseSystemConf, BaseSystemEvolve, BaseSystemCalc, BaseSystemPlot, BaseSystemGet):
    """
    The BaseSystem class is the base class for all systems in ComFiT.
    """
    
    def fft(self, field):
        """
        Perform a fast Fourier transform on a field.

        Parameters
        ----------
        field : np.ndarray
            The field to transform
            
        Returns
        -------
        np.ndarray
            The transformed field

        Note
        ----
        The field is assumed to represent a field in the dimensions of the system.
        Thus, if the system is 2D and the field is 3D, it is assumed that 
        the provided field is a collection of 2D fields which will be transformed
        individually.
        """
        
        return sp.fft.fftn(field, axes=range(-self.dim, 0))

    def ifft(self, field):
        """Perform an inverse fast Fourier transform on a field.

        Parameters
        ----------
        field : np.ndarray
            The field to transform
            
        Returns
        -------
        np.ndarray
            The transformed field

        Note
        ----
        The field is assumed to represent a field in the dimensions of the system.
        Thus, if the system is 2D and the field is 3D, it is assumed that 
        the provided field is a collection of 2D fields which will be transformed
        individually.
        """
        
        return sp.fft.ifftn(field, axes=range(-self.dim, 0))