###############################################################################
# demust - A Python toolkit to visualize and analyze deep mutational scanning #
#          data of proteins.                                                  #
# Authors: Mustafa Tekpinar                                                   #
# Copyright Mustafa Tekpinar 2022                                             #
#                                                                             #
# This file is part of demust.                                                #
#                                                                             #
# demust is free software: you can redistribute it and/or modify              #
# it under the terms of the GNU Lesser General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# demust is distributed in the hope that it will be useful,                   #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
# GNU LESSER General Public License for more details.                         #
#                                                                             #
# You should have received a copy of the GNU Lesser General Public License    #
# along with demust.  If not, see <https://www.gnu.org/licenses/>.            #
###############################################################################

from setuptools import setup, find_packages

from demust import __version__ as cp_vers

setup(name='demust',
      version=cp_vers,
      description="A Python toolkit to visualize and analyze deep mutational scanning data of proteins.",
      long_description=open('README.md').read(),
      long_description_content_type="text/markdown",
      author="Mustafa Tekpinar",
      author_email="tekpinar@buffalo.edu",
      url="https://github.com/tekpinar/demust",
      download_url="https://github.com/tekpinar/demust",
      license="LGPL",
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Environment :: Console',
          'Operating System :: POSIX',
          'Operating System :: Microsoft :: Windows',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Chemistry'
          ],
      python_requires='>=3.6',
      install_requires=[i for i in [l.strip() for l in open("requirements.txt").read().split('\n')] if i],
      # zip_safe=False,
      packages=[p for p in find_packages() if p != 'tests'],
      # file where some variables must be fixed by install
      entry_points={
          'console_scripts': [
              'demust=demust.scripts.demust:main'
          ]
      }
      )
