from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = Extension('solar_system_ext',
                sources=['solar_system_derivative.pyx', 'solar_system_data_cy.pxi'])

setup(name='spectral_norm',
      ext_modules=cythonize('spectral_norm.pyx'))
