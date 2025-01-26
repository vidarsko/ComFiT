import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import time

import plotly.graph_objects as go

bs = cf.BaseSystem(3, xRes=21, yRes=21, zRes=21)
bs.plot_lib = 'plotly'

vector_field = ([bs.x, bs.y, bs.z])

bs.plot_vector_field(vector_field, spacing=1).show()