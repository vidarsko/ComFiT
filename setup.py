from setuptools import setup

setup(
    name='comfit',
    version='1.8.1',
    packages=['comfit'],
    package_data={'comfit':['core/*',
                            'tool/*', 
                            'plot/*',
                            'quantum_mechanics/*', 
                            'bose_einstein_condensate/*',
                            'nematic_liquid_crystal/*',
                            'phase_field_crystal/*' ]},
    author='Vidar Skogvoll and Jonas Rønning',
    install_requires=['numpy>=1.21, <2.0.0',
                      'scipy',
                      'matplotlib',
                      'plotly',
                      'kaleido<=0.2.1',
                      'scikit-image>=0.18.0,<0.23',
                      'moviepy',
                      'pillow'],
)