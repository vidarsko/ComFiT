from comfit.core.base_system_init import BaseSystemInit
from comfit.core.base_system_conf import BaseSystemConf
from comfit.core.base_system_evolve import BaseSystemEvolve
from comfit.core.base_system_calc import BaseSystemCalc
from comfit.core.base_system_plot import BaseSystemPlot
from comfit.core.base_system_get import BaseSystemGet

import scipy as sp

class BaseSystem(BaseSystemInit, BaseSystemConf, BaseSystemEvolve, BaseSystemCalc, BaseSystemPlot, BaseSystemGet):
    """
    The BaseSystem class is the base class for all systems in ComFiT
    """
    
    def fft(self, field):
        """
        Perform a fast Fourier transform on a field. 

        Args: 
            field (array): The field to transform
        Returns:
            array: The transformed field
        """
        
        return sp.fft.fftn(field, axes=range(-self.dim, 0))

    def ifft(self, field):
        """
        Perform an inverse fast Fourier transform on a field.

        Args:
            field (array): The field to transform
        
        Return:
            array: The transformed field
        """
        
        return sp.fft.ifftn(field, axes=range(-self.dim, 0))