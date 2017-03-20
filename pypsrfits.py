#!/usr/bin/env python

# Copyright (C) 2017 by Paul Demorest and Maciej Serylak
# Licensed under the Academic Free License version 3.0
# This program comes with ABSOLUTELY NO WARRANTY.
# You are free to modify and redistribute this code as long
# as you do not remove the above attribution and reasonably
# inform recipients that you have modified the original work.

import numpy
import astropy.io.fits as pyfits

class PSRFITS:
  def __init__(self, fileName = None):
    self.fits = None
    if fileName != None:
      self.open(fileName)

  def open(self, fileName):
    """Open the specified PSRFITS file. A pyfits object is
    created and stored as self.fits. For convenience, the
    primary HDU table header is stored as self.hdr, and
    the subint HDU table header as self.subintHdr."""
    self.fits = pyfits.open(fileName, memmap = True, mode = "readonly", save_backup = False)
    self.hdr = self.fits[0].header
    self.subintHdr = self.fits["SUBINT"].header

  def getFreqs(self, row = 0):
    """Return the frequency array from the specified subint."""
    return self.fits["SUBINT"].data["DAT_FREQ"][row]

  def getData(self, startRow = 0, endRow = None, downSamp = 1, freqDownSamp = 1, applyScales = True, getFT = False, squeeze = False):
    """Read the data from the specified rows and return it as a
    single array with [time, poln, chan] dimensions.
    Options:
      startRow: first subint read (0-based index)
      endRow: final subint to read. None implies endRow = startRow.
              Negative values imply offset from the end, i.e.
              getData(0,-1) would read the entire file. (Don't forget
              that PSRFITS files are often huge so this might be a bad idea).
      downSamp: downsample the data in time as they are being read in.
                The downsample factor should evenly divide the number
                of spectra per row. downSamp = 0 means integrate each
                row completely.
      freqDownSamp: downsample the data in frequency as they are being read in.
                    The downsample factor should evenly divide the number
                    of channels.
      applyScales: set to False to avoid applying the scale/offset
                   data stored in the file.
      getFT: if True return time and frequency arrays as well.
      squeeze: if True, "squeeze" the data array (remove len-1
               dimensions).
    Notes: Only 8, 16, and 32 bit data are currently understood."""

    # Check if the mode of observations is correct.
    if self.hdr["OBS_MODE"].strip() != "SEARCH":
      raise RuntimeError("getData() only works on SEARCH-mode PSRFITS")

    # Get header parameters.
    nSampBlk = self.subintHdr["NSBLK"]
    nPol = self.subintHdr["NPOL"]
    nChan = self.subintHdr["NCHAN"]
    nBit = self.subintHdr["NBITS"]
    tBin = self.subintHdr["TBIN"]
    polType = self.subintHdr["POL_TYPE"]
    nRowsFile = self.subintHdr["NAXIS2"]
    #print "nSampBlk:", nSampBlk
    #print "nPol:", nPol
    #print "nChan:", nChan
    #print "nBit:", nBit
    #print "tBin:", tBin
    #print "polType:", polType
    #print "nRowsFile:", nRowsFile

    # Check if down sampling in time and frequency is enabled.
    if downSamp == 0:
      downSamp = 1
    if downSamp > nSampBlk:
      downSamp = nSampBlk
    if freqDownSamp == 0:
      freqDownSamp = 1
    if freqDownSamp > nChan:
      freqDownSamp = nChan
    if endRow == None:
      endRow = startRow
    if endRow < 0:
      endRow = nRowsFile + endRow
    if nSampBlk % downSamp > 0:
      raise RuntimeError("downSamp does not evenly divide NSBLK (%d)." % nSampBlk)
    if nChan % freqDownSamp > 0:
      raise RuntimeError("freqDownSamp does not evenly divide NCHAN (%d)." % nChan)

    nRowsTotal = endRow - startRow + 1
    downNSampBlk = nSampBlk / downSamp
    downNChan = nChan / freqDownSamp
    downTBin = tBin * downSamp
    #print "nRowsTotal", nRowsTotal
    #print "downNSampBlk", downNSampBlk
    #print "downNChan", downNChan
    #print "downTBin", downTBin

    # Data types of the signed and unsigned.
    if nBit == 8:
      signedType = numpy.int8
      unsignedType = numpy.uint8
    elif nBit == 16:
      signedType = numpy.int16
      unsignedType = numpy.uint16
    elif nBit == 32:
      signedType = numpy.float32
      unsignedType = numpy.float32
    else:
      raise RuntimeError("Not handled number of NBITS (%d)." % nBit)

    # Allocate arrays.
    rowSpectrumResult = numpy.zeros(nChan, dtype = numpy.float32)
    finalSpectrum = numpy.zeros((nRowsTotal * downNSampBlk, nPol, downNChan), dtype = numpy.float32)
    if getFT:
      finalFreqs = numpy.zeros(downNChan)
      finalTiming = numpy.zeros(nRowsTotal * downNSampBlk)
    polSign = 1
    if "AABB" in polType:
      polSign = 2

    # Iterate over rows (blocks).
    for iRow in range(nRowsTotal):
      if applyScales:
        offsets = self.fits["SUBINT"].data["DAT_OFFS"][iRow + startRow]
        scales = self.fits["SUBINT"].data["DAT_SCL"][iRow + startRow]
        weights = self.fits["SUBINT"].data["DAT_WTS"][iRow + startRow]
        offsets = offsets.reshape((nPol, nChan))
        scales = scales.reshape((nPol, nChan))
        #weights = numpy.concatenate((weights, weights,weights, weights))
        #weights = weights.reshape((nPol, nChan))
      if getFT:
        rowTiming = self.fits["SUBINT"].data["OFFS_SUB"][iRow + startRow] - (self.fits["SUBINT"].data["TSUBINT"][iRow + startRow] / 2.0)
        freqs_row = self.fits["SUBINT"].data["DAT_FREQ"][iRow + startRow]
      rowData = self.fits["SUBINT"].data["DATA"][iRow + startRow]
      #print "rowData", rowData.shape
      #print "offets", offsets.shape
      #print "scales", scales.shape
      #print "weights", weights.shape

      # Decode rowDAta for 16 bit data type file.
      if (nBit == 16):
        rowData = numpy.fromstring(rowData.tostring(), dtype = numpy.int16)
        rowData = rowData.reshape((nSampBlk, nPol, nChan))

      # Iterate over samples and polarisations.
      for iSamp in range(downNSampBlk):
        for iPol in range(nPol):
          if iPol < polSign:
            dataType = unsignedType
          else:
            dataType = signedType
          rowSpectrumResult = rowData[iSamp * downSamp:(iSamp + 1) * downSamp, iPol, :].astype(dataType).mean(0)
          if getFT:
            finalTiming[iRow * downNSampBlk + iSamp] = rowTiming + (iSamp + 0.5) * downTBin
          if applyScales:
            rowSpectrumResult *= scales[iPol, :]
            rowSpectrumResult += offsets[iPol, :]
            rowSpectrumResult *= weights[iPol, :]
          if freqDownSamp == 1:
            finalSpectrum[iRow * downNSampBlk + iSamp, iPol, :] = rowSpectrumResult
            if getFT:
              finalFreqs[:] = freqs_row[:]
          else:
            finalSpectrum[iRow * downNSampBlk + iSamp, iPol, :] = rowSpectrumResult.reshape((-1, freqDownSamp)).mean(1)
            if getFT:
              finalFreqs[:] = freqs_row.reshape((-1, freqDownSamp)).mean(1)

    if squeeze:
      finalSpectrum = finalSpectrum.squeeze()
    if getFT:
      return (finalSpectrum, finalTiming, finalFreqs)
    else:
      return finalSpectrum
