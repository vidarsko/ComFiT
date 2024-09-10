import numpy as np

def tool_complete_field(self, field: np.ndarray) -> np.ndarray:
        """Extends the field in case not a complete array is given. 

        For instance, if a field in 3 dimensions is calculated using only x, then the field is extended to the full 3D array.

        Args:
            field: The field to be extended

        Returns:
            The extended field (np.ndarray)
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