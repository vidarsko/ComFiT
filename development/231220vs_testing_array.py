import numpy as np

# Example field array
field = np.random.rand(55, 55, 55)

# Example index arrays
R0 = np.random.randint(0, 55, size=(216, 216))
R1 = np.random.randint(0, 55, size=(216, 216))
R2 = np.random.randint(0, 55, size=(216, 216))

# Extract values
result = field[R0, R1, R2]