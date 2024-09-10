import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

nem = cf.NematicLiquidCrystal(2,xRes=64,yRes=64,dx=1,dy=1,dt=0.1)



nem.conf_initial_condition_disordered(noise_strength=4)
nem.evolve_nematic_no_flow(10)

for i in range(40):
    nem.evolve_nematic(40)
    test_1 =  nem.Q_f[0][1]
    test_2 = nem.Q_f[0,1,:,:]
    tester1 = np.abs(test_1 - test_2)

    print(tester1.sum())
