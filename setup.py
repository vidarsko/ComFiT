from setuptools import setup

setup(
    name='comfit',
    version='1.1.0',
    packages=['comfit'],
    package_data={'comfit':['core/*','models/*','tools/*']},
    author='Vidar Skogvoll and Jonas Rønning',
    install_requires=['numpy>=1.22.0',
                      'scipy',
                      'scikit-image',
                      'matplotlib',
                      'moviepy==1.0.3',
                      'imageio',
                      'apptools==5.2.1',
                      'envisage==7.0.3',
                      'packaging==23.1',
                      'traits==6.4.3',
                      'pyface==8.0.0',
                      'pygments==2.16.1',
                      'traitsui==8.0.0',
                      'vtk==9.2.6',
                      'mayavi==4.8.1'],
)