import numpy as np
import matplotlib.pyplot as plt

Q = np.zeros((2,2,100,100))

Q[0,0,:,:] =1
v = np.random.randint(4,size=(2,100,100))
alpha = np.random.randint(4,size=(100,100))

beta = np.einsum('iklm,klm-> ilm',Q,v)

print()