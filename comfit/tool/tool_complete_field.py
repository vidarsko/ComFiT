from typing import Union, Tuple, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from comfit.core.base_system import BaseSystem

def tool_complete_field(
        self : 'BaseSystem', 
        field : np.ndarray
        ) -> np.ndarray:
    """Extends a partial field array to a complete array by tiling.

    Parameters
    ----------
    self : 'BaseSystem'
        A BaseSystem (or derived) instance.
    field : np.ndarray
        Input field array to be extended. Can be 1D, 2D, or 3D array with partial dimensions.

    Returns
    -------
    np.ndarray
        Extended field array with complete dimensions matching self.xRes, self.yRes, and self.zRes.

    Notes
    -----
    If input field is incomplete (has dimensions of 1 in some axes), it will be tiled
    to match the full resolution specified by self.xRes, self.yRes, and self.zRes.
    """

    # 2 dimensional fields
    if field.shape == (self.xRes,1):
        field = np.tile(field,(1,self.yRes))

    elif field.shape == (1,self.yRes):
        field = np.tile(field,(self.xRes,1))

    # 3 dimensional fields
    elif field.shape == (self.xRes,1,1):
        field = np.tile(field,(1,self.yRes,1))
        field = np.tile(field,(1,1,self.zRes))

    elif field.shape == (1,self.yRes,1):
        field = np.tile(field,(self.xRes,1,1))
        field = np.tile(field,(1,1,self.zRes))

    elif field.shape == (1,1,self.zRes):
        field = np.tile(field,(self.xRes,1,1))
        field = np.tile(field,(1,self.yRes,1))

    elif field.shape == (self.xRes,self.yRes,1):
        field = np.tile(field,(1,1,self.zRes))

    elif field.shape == (self.xRes,1,self.zRes):
        field = np.tile(field,(1,self.yRes,1))

    elif field.shape == (1,self.yRes,self.zRes):
        field = np.tile(field,(self.xRes,1,1))

    return field