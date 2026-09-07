from setuptools import setup

setup(
    name='comfit',
    version='1.9.6',
    packages=['comfit'],
    package_data={'comfit':['core/*',
                            'tool/*', 
                            'plot/*',
                            'quantum_mechanics/*', 
                            'bose_einstein_condensate/*',
                            'nematic_liquid_crystal/*',
                            'phase_field_crystal/*' ]},
    author='Vidar Skogvoll and Jonas Rønning',
    install_requires=['numpy',
                      'scipy',
                      'matplotlib',
                      'plotly>=6.1.1',
                      'kaleido>=1.0.0',
                      'scikit-image>=0.25.2',
                      'moviepy',
                      'pillow'],
)