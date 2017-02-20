#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name="pypsrfits",
      description="Simple code to read search-mode PSRFITS data into numpy array.",
      author="Paul Demorest and Maciej Serylak",
      author_email="mserylak@ska.ac.za",
      packages=find_packages(),
      url='https://github.com/mserylak/pypsrfits',
      license="Academic Free License 3.0",
      classifiers=[
          "Development Status :: 1 - Beta",
          "Intended Audience :: Developers",
          "License :: OSI Approved :: AFL-3.0 License",
          "Operating System :: OS Independent",
          "Programming Language :: Python",
          "Programming Language :: Python :: 2",
          "Programming Language :: Python :: 2.6",
          "Programming Language :: Python :: 2.7",
          "Topic :: Software Development :: Libraries :: Python Modules",
          "Topic :: Scientific/Engineering :: Astronomy"],
      platforms=["OS Independent"],
      keywords="psrfits",
      zip_safe=False,
      install_requires=['numpy', 'pyfits'])
