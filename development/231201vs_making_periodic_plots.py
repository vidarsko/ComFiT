import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

s = cf.BaseSystem(1)

r2 = (s.x-s.xmid)**2


R2 = (0.4*s.xmax)**2


ax = plt.gcf().add_subplot(111)

w2 = (0.2*s.xmax)**2
y = 1/2*(1+np.tanh((r2-R2)/w2))
ax.plot(s.x,y,label='$w=0.2 x_{max}$')

w2 = (0.1*s.xmax)**2
y = 1/2*(1+np.tanh((r2-R2)/w2))
ax.plot(s.x,y,label='$w=0.1 x_{max}$')

ax.legend()
ax.grid(True)
ax.set_xlabel('$x/a_0$')
ax.set_ylabel('$F$')

plt.gcf().set_size_inches(2.5,2.5)
plt.tight_layout()
plt.savefig(r'C:\Users\vidar\UiO Dropbox\Vidar Skogvoll\Apps\Overleaf\230531 - Code documentation\Figures\NumericalImplementationOfAngleFields\FilterFunction.png',dpi=300)