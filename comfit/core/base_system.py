from comfit.core.base_system_init import BaseSystemInit
from comfit.core.base_system_conf import BaseSystemConf
from comfit.core.base_system_evolve import BaseSystemEvolve
from comfit.core.base_system_calc import BaseSystemCalc
from comfit.core.base_system_plot import BaseSystemPlot
from comfit.core.base_system_get import BaseSystemGet


class BaseSystem(BaseSystemInit, BaseSystemConf, BaseSystemEvolve, BaseSystemCalc, BaseSystemPlot, BaseSystemGet):
    """
    The BaseSystem class is the base class for all systems in ComFiT.
    """
    pass