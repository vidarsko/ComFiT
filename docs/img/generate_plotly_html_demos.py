"""Regenerate the interactive Plotly HTML embeds used in the documentation.

Run from the repository root with `python docs/img/generate_plotly_html_demos.py`.
Each embed is a self-contained HTML file (Plotly.js loaded from a CDN) that gets
included in a docs page via an <iframe>. The code in each function below is kept
in sync with the corresponding `=== "plotly"` tab in docs/Plotting.md.
"""

import numpy as np
import comfit as cf

OUT_DIR = 'docs/img'


def _save(fig, filename, height=420):
    # Fixed light background so the embed stays legible regardless of the
    # docs site's light/dark theme (the interactive plot isn't theme-aware).
    fig.update_layout(
        height=height,
        margin=dict(l=10, r=10, t=30, b=10),
        paper_bgcolor='white',
        plot_bgcolor='white',
    )
    fig.write_html(f'{OUT_DIR}/{filename}', include_plotlyjs='cdn', full_html=True)
    print(f'Wrote {OUT_DIR}/{filename}')


def generate_plot_field_demo():
    bs1 = cf.BaseSystem(1, xRes=31)
    field1 = bs1.x**2

    bs2 = cf.BaseSystem(2, xRes=31, yRes=31)
    field2 = bs2.x**2 + bs2.y**2

    bs3 = cf.BaseSystem(3, xRes=31, yRes=31, zRes=31)
    field3 = bs3.x**2 + bs3.y**2 + bs3.z**2

    fig, axs = bs1.plot_subplots(1, 3)

    bs1.plot_field(field1, fig=fig, ax=axs[0])
    bs2.plot_field(field2, fig=fig, ax=axs[1])
    bs3.plot_field(field3, fig=fig, ax=axs[2])

    _save(fig, 'plotting_plot_field_demo_interactive.html')


def generate_plot_field_in_plane_demo():
    bs = cf.BaseSystem(3, xRes=31, yRes=31, zRes=31)
    field = bs.x**2 + bs.y**2 + bs.z**2

    fig, axs = bs.plot_subplots(1, 2)

    bs.plot_field_in_plane(field, fig=fig, ax=axs[0])
    bs.plot_field_in_plane(field, fig=fig, ax=axs[1], normal_vector=[1, 1, 0], position=[10, 10, 10])

    _save(fig, 'plotting_plot_field_in_plane_demo_interactive.html')


def generate_plot_complex_field_demo():
    bs1 = cf.BaseSystem(1, xRes=31)
    field1 = bs1.x**2 * np.exp(1j * bs1.x / 3)

    bs2 = cf.BaseSystem(2, xRes=31, yRes=31)
    field2 = (bs2.x**2 + bs2.y**2) * np.exp(1j * bs2.x / 3)

    # Lower resolution than the matplotlib example: the 3D phase_angle/phase_blob
    # methods each add several marching-cubes surfaces, which makes the embedded
    # file grow quickly with resolution.
    bs3 = cf.BaseSystem(3, xRes=15, yRes=15, zRes=15)
    field3 = (bs3.x**2 + bs3.y**2 + bs3.z**2) * np.exp(1j * bs3.x / 3)

    fig, axs = bs1.plot_subplots(2, 3)

    bs1.plot_complex_field(field1, fig=fig, ax=axs[0][0])
    bs2.plot_complex_field(field2, fig=fig, ax=axs[0][1], plot_method='phase_angle')
    bs2.plot_complex_field(field2, fig=fig, ax=axs[0][2], plot_method='3Dsurface')
    bs3.plot_complex_field(field3, fig=fig, ax=axs[1][1], plot_method='phase_angle')
    bs3.plot_complex_field(field3, fig=fig, ax=axs[1][2], plot_method='phase_blob')

    _save(fig, 'plotting_plot_complex_field_demo_interactive.html', height=560)


def generate_plot_complex_field_in_plane_demo():
    bs = cf.BaseSystem(3, xRes=31, yRes=31, zRes=31)
    complex_field = (bs.x**2 + bs.y**2 + bs.z**2) * np.exp(1j * bs.y / 3)

    fig, axs = bs.plot_subplots(1, 2)

    bs.plot_complex_field_in_plane(complex_field, fig=fig, ax=axs[0])
    bs.plot_complex_field_in_plane(complex_field, fig=fig, ax=axs[1], normal_vector=[0, 0, 1], position=[10, 10, 10])

    _save(fig, 'plotting_plot_complex_field_in_plane_demo_interactive.html')


def generate_plot_angle_field_demo():
    bs1 = cf.BaseSystem(1, xRes=31)
    angle_field1 = np.mod(bs1.x / 5, 2 * np.pi) - np.pi

    bs2 = cf.BaseSystem(2, xRes=31, yRes=31)
    angle_field2 = np.mod((bs2.x + 2 * bs2.y) / 5, 2 * np.pi) - np.pi

    # Lower resolution than the matplotlib example, for the same reason as
    # generate_plot_complex_field_demo above (dim=3 uses the same isosurface method).
    bs3 = cf.BaseSystem(3, xRes=15, yRes=15, zRes=15)
    angle_field3 = np.mod((bs3.x + 2 * bs3.y + 3 * bs3.z) / 5, 2 * np.pi) - np.pi

    fig, axs = bs1.plot_subplots(1, 3)

    bs1.plot_angle_field(angle_field1, fig=fig, ax=axs[0])
    bs2.plot_angle_field(angle_field2, fig=fig, ax=axs[1])
    bs3.plot_angle_field(angle_field3, fig=fig, ax=axs[2])

    _save(fig, 'plotting_plot_angle_field_demo_interactive.html')


def generate_plot_angle_field_in_plane_demo():
    bs = cf.BaseSystem(3, xRes=31, yRes=31, zRes=31)
    angle_field = np.mod((bs.x + 2 * bs.y + 3 * bs.z) / 5, 2 * np.pi) - np.pi

    fig, axs = bs.plot_subplots(1, 2)

    bs.plot_angle_field_in_plane(angle_field, fig=fig, ax=axs[0])
    bs.plot_angle_field_in_plane(angle_field, fig=fig, ax=axs[1], normal_vector=[0, 0, 1], position=[10, 10, 10])

    _save(fig, 'plotting_plot_angle_field_in_plane_demo_interactive.html')


def generate_plot_vector_field_demo():
    bs1 = cf.BaseSystem(1, xRes=31)
    bs2 = cf.BaseSystem(2, xRes=31, yRes=31)
    bs3 = cf.BaseSystem(3, xRes=11, yRes=11, zRes=11)

    fig, axs = bs1.plot_subplots(3, 3)

    vector_field = np.array([bs1.x * np.cos(bs1.x / 5)])
    bs1.plot_vector_field(vector_field, fig=fig, ax=axs[0][0], spacing=1)

    vector_field = np.array([bs1.x * np.cos(bs1.x / 5), bs1.x * np.sin(bs1.x / 5)])
    bs1.plot_vector_field(vector_field, fig=fig, ax=axs[0][1], spacing=2)

    vector_field = np.array([bs1.x * np.cos(bs1.x / 5), bs1.x * np.sin(bs1.x / 5), bs1.x * np.cos(bs1.x / 5)])
    bs1.plot_vector_field(vector_field, fig=fig, ax=axs[0][2], spacing=3)

    vector_field = np.array([bs2.x * np.cos(bs2.y / 5)])
    bs2.plot_vector_field(vector_field, fig=fig, ax=axs[1][0], spacing=3)

    vector_field = np.array([bs2.x * np.cos(bs2.y / 5), bs2.y * np.sin(bs2.x / 5)])
    bs2.plot_vector_field(vector_field, fig=fig, ax=axs[1][1], spacing=5)

    vector_field = np.array([bs2.x * np.cos(bs2.y / 5), bs2.y * np.sin(bs2.x / 5), bs2.x * np.cos(bs2.y / 5)])
    bs2.plot_vector_field(vector_field, fig=fig, ax=axs[1][2], spacing=3)

    vector_field = np.array([bs3.z + bs3.x * np.cos(bs3.y / 5)])
    bs3.plot_vector_field(vector_field, fig=fig, ax=axs[2][0], spacing=3)

    vector_field = np.array([bs3.z + bs3.x * np.cos(bs3.y / 5), bs3.z + bs3.y * np.sin(bs3.x / 5)])
    bs3.plot_vector_field(vector_field, fig=fig, ax=axs[2][1], spacing=5)

    vector_field = np.array([
        bs3.z + bs3.x * np.cos(bs3.y / 5),
        bs3.z + bs3.y * np.sin(bs3.x / 5),
        -bs3.z + bs3.x * np.cos(bs3.y / 5),
    ])
    bs3.plot_vector_field(vector_field, fig=fig, ax=axs[2][2], spacing=3)

    _save(fig, 'plotting_plot_vector_field_demo_interactive.html', height=700)


def generate_plot_vector_field_in_plane_demo():
    bs = cf.BaseSystem(3, xRes=11, yRes=11, zRes=11)

    fig, axs = bs.plot_subplots(1, 2)

    vector_field = np.array([
        bs.z + bs.x * np.cos(bs.y / 5),
        bs.z + bs.y * np.sin(bs.x / 5),
        -bs.z + bs.x * np.cos(bs.y / 5),
    ])
    bs.plot_vector_field_in_plane(vector_field, fig=fig, ax=axs[0])

    vector_field2 = np.array([bs.z + bs.x * np.cos(bs.y / 5), bs.z + bs.y * np.sin(bs.x / 5)])
    bs.plot_vector_field_in_plane(vector_field2, fig=fig, ax=axs[1], normal_vector=[0, 1, 1], position=[2, 3, 3])

    _save(fig, 'plotting_plot_vector_field_in_plane_demo_interactive.html')


if __name__ == '__main__':
    generate_plot_field_demo()
    generate_plot_field_in_plane_demo()
    generate_plot_complex_field_demo()
    generate_plot_complex_field_in_plane_demo()
    generate_plot_angle_field_demo()
    generate_plot_angle_field_in_plane_demo()
    generate_plot_vector_field_demo()
    generate_plot_vector_field_in_plane_demo()
