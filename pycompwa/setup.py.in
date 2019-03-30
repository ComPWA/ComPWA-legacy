#!/usr/bin/env python
from setuptools import setup, find_packages

setup(name='pycompwa',
      version='0.0.1',
      author='The ComPWA team',
      maintainer = "Peter Weidenkaff",
      maintainer_email = "weidenka@uni-mainz.de",
      url = "https://github.com/ComPWA/ComPWA",
      description='ComPWA: The common Partial Wave Analysis framework',
      long_description='',
      license = "GPLv3 or later",
      package_dir={'': '${CMAKE_CURRENT_BINARY_DIR}'},
      packages=find_packages(),
      package_data={
        # Include default particle list and precompiled pybind interface
        '': ['particle_list.xml', 'ui*'],
        # And include any *.msg files found in the 'hello' package, too:
        #  'hello': ['*.msg'],
      },
      zip_safe=False,
      test_suite = "tests",
      #  test_requires=['numpy>=1.14.5', 'pytest>=3.6.3',
      #                    'xmltodict>=0.11.0', 'scipy>=1.1.0',
      #                    'uproot>=3.2.5', 'matplotlib>=2.2.2'],
      install_requires=['numpy>=1.14.5', 'pytest>=3.6.3',
                        'xmltodict>=0.11.0', 'scipy>=1.1.0',
                        'uproot>=3.2.5', 'matplotlib>=2.2.2'],
      )
