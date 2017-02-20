# Introduction

The original code is authored by Paul Demorest and can be found in https://github.com/demorest/pypsrfits

pypsrfits is a very simple python module for reading search-mode PSRFITS data
into numpy array.

# Dependencies

The module requires the pyfits and numpy python modules.

# Example usage

Importing and loading data file into the pyfits object:

`import pypsrfits`

`f = pypsrfits.PSRFITS("my_file.fits")`

The full pyfits object for the file is available:

`f.fits`

The main header and SUBINT header are also accessible:

`f.hdr`

`f.subintHdr`

To read all of the data from row 13:

`d = f.getData(13)`

Read all data in entire file, downsampling in time by
a factor of 256

`d = f.getData(0, -1, downSamp = 256)`
