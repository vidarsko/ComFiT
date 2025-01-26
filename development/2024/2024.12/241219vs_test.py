import numpy as np
import matplotlib.pyplot as plt

def f(t, y):
    # Example ODE: dy/dt = -y
    return -y

def forward_euler(f, y0, t0, tf, dt):
    t = np.arange(t0, tf, dt)
    y = np.zeros(len(t))
    y[0] = y0
    
    for i in range(1, len(t)):
        y[i] = y[i-1] + dt * f(t[i-1], y[i-1])
    
    return t, y

# Initial conditions
y0 = 1.0
t0 = 0.0
tf = 5.0
dt = 0.1

# Solve the ODE
t, y = forward_euler(f, y0, t0, tf, dt)

# Plot the solution
plt.plot(t, y, label='Forward Euler')
plt.xlabel('Time t')
plt.ylabel('y(t)')
plt.title('Solution of dy/dt = -y using Forward Euler Method')
plt.legend()
plt.grid(True)
plt.show()