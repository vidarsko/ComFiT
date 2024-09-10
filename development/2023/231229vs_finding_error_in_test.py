import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint

pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(nx=13, ny=13, nz=13)
# Set the plane equation
a = 0.3
b = 0.7
c = 0.1

for dislocation_type in range(2,3):
    # Insert a dislocation ring
    eta = pfc.calc_amplitudes_with_dislocation_ring(
    normal_vector=[a,b,c],
    position=pfc.rmid,
    radius=pfc.xmax/3,
    dislocation_type=dislocation_type
    )
    pfc.conf_PFC_from_amplitudes(eta)

    # Relax the system
    pfc.evolve_PFC(20)
    dislocation_nodes = pfc.calc_dislocation_nodes()


    #Check that all the nodes belong to the correct plane
    for dislocation_node in dislocation_nodes:


        # Check that the Burgers vector is correct
        print(dislocation_node['Burgers_vector'])
        print(pfc.a[dislocation_type-1])
        print((dislocation_node['Burgers_vector'] == pfc.a[dislocation_type-1]).all())

