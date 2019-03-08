#!/usr/bin/env python
from setuptools import setup, find_packages

setup(name='pycompwa',
      version='0.0.1',
      author='Hans Gans',
      description='A test project using pybind11 and CMake',
      long_description='',
      #  package_dir={ '': '${CMAKE_CURRENT_SOURCE_DIR}' },
      #  package_dir={ '': '${CMAKE_CURRENT_SOURCE_DIR}' },
      #  packages=['pycompwa'],
      packages=find_packages(),
      package_data={
        '': ['particle_list.xml'],
        # And include any *.msg files found in the 'hello' package, too:
        #  'hello': ['*.msg'],
      },
      zip_safe=False,
      install_requires=['numpy>=1.14.5', 'pytest>=3.6.3',
                        'xmltodict>=0.11.0', 'scipy>=1.1.0',
                        'uproot>=3.2.5', 'matplotlib>=2.2.2'],
)

