####################################################################################################
# # Copyright (C) 2024-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent

# Read requirements text file to get required packages
with (here / 'requirements.txt').open() as reqsfile:
    reqs = [
        line.strip()
        for line in reqsfile
        if line.strip() and not line.startswith('#')
    ]

# Setup package with params
setup(name='pyscev',
      version='0.2.0',
      description='pySCEv: pySCE single cell entropy scoring plus angular RNA velocity entropy',
      url='https://github.com/oliviahelens/pySCEv',
      author='Daniel Charytonowicz (pySCE); pySCEv fork contributors',
      author_email='daniel.charytonowicz@icahn.mssm.edu',
      license='GNU GPLv3',
      install_requires=reqs,
      extras_require={
          'velocity': ['scvelo>=0.2.4'],
          'test': ['pytest'],
      },
      packages=find_packages(include=['pysce', 'pysce.*']),
      package_data={'pysce': ['data/ppi_scent_pp_sparse.h5ad']},
      include_package_data = True,
      zip_safe=False)
