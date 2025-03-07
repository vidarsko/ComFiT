import comfit as cf
import matplotlib.pyplot as plt
import time




cfi = cf.BaseSystem(2, xRes=64, yRes=64)
field = cfi.x + cfi.y

times = {}

for plot_lib in ['plotly', 'matplotlib']:
    start_time = time.time()
    
    for n in range(100):  
        fig, ax = cfi.plot_field(field, plot_lib=plot_lib)
        cfi.plot_save(fig, n)  

    cf.tool_make_animation_gif(n)
    
    times[plot_lib] = time.time() - start_time

print(f"Plotly time: {times['plotly']:.2f}s")
print(f"Matplotlib time: {times['matplotlib']:.2f}s")
print(f"Ratio (Plotly/Matplotlib): {times['plotly']/times['matplotlib']:.2f}")

