from setuptools import setup

setup(
    name='comfit',
    version='1.3.0',
    packages=['comfit'],
    package_data={'comfit':['core/*','models/*','tools/*']},
    author='Vidar Skogvoll and Jonas Rønning',
    install_requires=['numpy',
                      'scipy',
                      'matplotlib',
                      'scikit-image',
                      'moviepy',
                      'imageio'],
)