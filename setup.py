from setuptools import setup

setup(
    name='comfit',
    version='1.1.0',
    packages=['comfit'],
    package_data={'comfit':['core/*','models/*','tools/*']},
    author='Vidar Skogvoll and Jonas Rønning',
    install_requires=['numpy',
                      'scipy',
                      'scikit-image',
                      'matplotlib',
                      'moviepy==1.0.3',
                      'imageio',
                      'vtk==9.2.6',
                      'mayavi==4.8.1',
                      'PyQt5==5.12.2'],
)