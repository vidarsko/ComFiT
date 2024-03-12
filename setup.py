from setuptools import setup

setup(
    name='comfit',
    version='1.3.1',
    packages=['comfit'],
    package_data={'comfit':['core/*','models/*','tools/*']},
    author='Vidar Skogvoll and Jonas Rønning',
    install_requires=['numpy>=1.21, <2.0.0',
                      'scipy',
                      'matplotlib',
                      'scikit-image',
                      'moviepy',
                      'imageio'],
)