# from distutils.core import setup
from setuptools import setup, find_packages
from PePr import __version__
setup(name="PePr",
      version=__version__, # change the version info in the PePr.__init__ file
      author="Yanxiao Zhang",
      author_email="troublezhang@gmail.com",
      url="https://github.com/shawnzhangyx/PePr/",
      license="GNU GPL v3",
      description="Peak-calling and Prioritization pipeline for replicated ChIP-Seq data",
      long_description="Peak-calling and Prioritization pipeline for replicated ChIP-Seq data",
      platforms = ['any'],
      classifiers=[
          "Development Status :: 4 - Beta",
          "Intended Audience :: Developers",
          "Intended Audience :: Science/Research",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
          "Programming Language :: Python",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
          ],
      packages=find_packages(),
      package_data={"PePr": ['data/*.bed']},
      install_requires=[
          'numpy>=1.6.0',
          'scipy>=0.14.0',
          'pysam',
          'sharedmem',
          ],
      entry_points={
          'console_scripts': [
              'PePr=PePr.PePr:argless_main',
              ]
          }
      )
