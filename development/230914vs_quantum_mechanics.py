import comfit as cf
import matplotlib.pyplot as plt

qm = cf.QM(1)
qm.set_initial_condition_gaussian()
qm.plot()
plt.show()
