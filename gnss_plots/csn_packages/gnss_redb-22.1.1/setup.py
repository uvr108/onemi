# -*- coding: utf-8 -*-
"""
*gnss_redb* is a package that allows to get data from a RethinkDB databases
with GNSS data. Designed for internal use in the CSN.

Copyright (c) 2020 Centro Sismológico Nacional, Universidad de Chile

This file is part of "gnss_redb".

"gnss_redb" is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

"gnss_redb" is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with "gnss_redb". If not, see <https://www.gnu.org/licenses/>.
"""
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
from gnss_redb import __version__
here = path.abspath(path.dirname(__file__))
package_name = 'gnss_redb'

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name=package_name,
    version=__version__,  # YYYYY.MM = year and month;
    description='Data Base convenience classes and methods (CSN internal use)',
    long_description=long_description,  # [1]
    # The project's main homepage.
    # url='https://github.com/pypa/sampleproject',
    # Author details
    author='Francisco del Campo R.',
    author_email='fdelcampo@csn.uchile.cl',
    license='GPLv3',  # [2]
    python_requires='>=3.5',
    install_requires=['six',
                      'rethinkdb',
                      'numpy',
                      'python-dateutil'],
    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    # packages=find_packages(exclude=['contrib', 'docs', 'examples']),
    packages=find_packages(),
    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_dir={package_name: package_name},
    package_data={package_name: ['../ref_coords/gnss_network'
                                 ]},
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        # 'Development Status :: 3 - Alpha',  # [3]
        # 'Development Status :: 4 - Beta',  # [3]
        'Development Status :: 5 - Production/Stable',  # [3]
        # Indicate who your project is intended for
        'Intended Audience :: Science/Research :: Developers',
        'Topic :: Scientific/Engineering',
        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 3 or both.
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8'
        ],
    # What does your project relate to?
    keywords='',
    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    # entry_points={
    #     'console_scripts': [
    #         'example=examples:main',
    #     ],
    # },
    )

# [1] on the PyPI the field "long_description" will be used
#     as the description of the package on the website.
# [2] Specifying a licence is important for Open Source software
# [3] Maturity of package.
