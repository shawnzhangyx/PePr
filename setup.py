#from distutils.core import setup
from setuptools import setup
setup(name="PePr",
      version="1.0.7",
      description="Peak-calling and Prioritization pipeline for replicated ChIP-Seq data",
      author="Yanxiao Zhang",
      author_email="troublezhang@gmail.com",
      url="https://code.google.com/p/pepr-chip-seq/",
      license="GNU GPL v3",
      classifiers=[
          "Development Status :: 4 - Beta", 
          "Intended Audience :: Developers",
          "Intended Audience :: Science/Research",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
          "Programming Language :: Python",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
          ],
      packages=["PePr"],
      package_data={"PePr":['data/*.bed']},
      install_requires=[
          'numpy==1.6.0',
          'scipy==0.14.0',
          'pysam',
          ]
      )



