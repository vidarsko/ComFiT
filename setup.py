from setuptools import setup

setup(
    name='comfit',
    version='1.0.0',
    packages=['comfit'],
    package_data={'comfit':['core/*','models/*']},
    author='Vidar Skogvoll and Jonas Rønning',
    install_requires=['numpy>=1.22.0','scikit-image']
)