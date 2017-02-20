#!/usr/bin/env python

from setuptools import setup

setup(name="pypsrfits",
      version="1.0",
      description="Simple code to read search-mode PSRFITS data into numpy array.",
      author="Paul Demorest and Maciej Serylak",
      author_email="mserylak@ska.ac.za",
      py_modules=["pypsrfits"],
      url='https://github.com/mserylak/pypsrfits',
      license="Academic Free License 3.0",
      classifiers=[
          "Development Status :: 4 - Beta",
          "Intended Audience :: Developers",
          "License :: OSI Approved :: Academic Free License (AFL)",
          "Operating System :: OS Independent",
          "Programming Language :: Python",
          "Programming Language :: Python :: 2.7",
          "Topic :: Software Development :: Libraries :: Python Modules",
          "Topic :: Scientific/Engineering :: Astronomy"],
      platforms=["OS Independent"],
      keywords="pulsar psrfits meerkat ska",
      zip_safe=False,
      install_requires=['numpy', 'pyfits'])
