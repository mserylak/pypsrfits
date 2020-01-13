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
    def __init__(self, filename=None):
        self.fits = None
        if filename is not None:
            self.open(filename)

    def open(self, filename):
        """Open the specified PSRFITS file. A pyfits object is
        created and stored as self.fits. For convenience, the
        primary HDU table header is stored as self.hdr, and
        the subint HDU table header as self.subintHdr."""
        self.fits = pyfits.open(filename, memmap=True, mode='readonly', save_backup=False)
        self.hdr = self.fits[0].header
        self.subint_hdr = self.fits['SUBINT'].header

    def get_freqs(self, row=0):
        """Return the frequency array from the specified subint."""
        return self.fits['SUBINT'].data['DAT_FREQ'][row]

    def get_data(self, start_row=0, end_row=None, downsamp=1, freq_downsamp=1, apply_scales=True, get_ft=False, squeeze=False):
        """Read the data from the specified rows and return it as a
        single array with [time, poln, chan] dimensions.
        Options:
            start_row: first subint read (0-based index)
            end_row: final subint to read. None implies endRow = startRow.
                     Negative values imply offset from the end, i.e.
                     getData(0,-1) would read the entire file. (Don't forget
                     that PSRFITS files are often huge so this might be a bad idea).
            down_samp: downsample the data in time as they are being read in.
                       The downsample factor should evenly divide the number
                       of spectra per row. downSamp = 0 means integrate each
                       row completely.
            freq_down_samp: downsample the data in frequency as they are being read in.
                            The downsample factor should evenly divide the number
                            of channels.
            apply_scales: set to False to avoid applying the scale/offset
                          data stored in the file.
            get_ft: if True return time and frequency arrays as well.
            squeeze: if True, "squeeze" the data array (remove len-1
                     dimensions).
         Notes: Only 8, 16, and 32 bit data are currently understood."""

        # Check if the mode of observations is correct.
        if self.hdr['OBS_MODE'].strip() != 'SEARCH':
            raise RuntimeError('get_data() only works on SEARCH-mode PSRFITS')

        # Get header parameters.
        nsamp_blk = self.subint_hdr['NSBLK']
        npol = self.subint_hdr['NPOL']
        nchan = self.subint_hdr['NCHAN']
        nbit = self.subint_hdr['NBITS']
        tbin = self.subint_hdr['TBIN']
        pol_type = self.subint_hdr['POL_TYPE']
        nrows_file = self.subint_hdr['NAXIS2']
        # print('nsamp_blk: {0}'.format(nsamp_blk))
        # print('npol: {0}'.format(npol))
        # print('nchan: {0}'.format(nchan))
        # print('nbit: {0}'.format(nbit))
        # print('tbin: {0}'.format(tbin))
        # print('pol_type: {0}'.format(pol_type))
        # print('nrows_file: {0}'.format(nrows_file))

        # Check if down sampling in time and frequency is enabled.
        if downsamp == 0:
            downsamp = 1
        if downsamp > nsamp_blk:
            downsamp = nsamp_blk
        if freq_downsamp == 0:
            freq_downsamp = 1
        if freq_downsamp > nchan:
            freq_downsamp = nchan
        if end_row is None:
            end_row = start_row
        if end_row < 0:
            end_row = nrows_file + end_row
        if nsamp_blk % downsamp > 0:
            raise RuntimeError('downsamp does not evenly divide NSBLK ({0}).'.format(nsamp_blk))
        if nchan % freq_downsamp > 0:
            raise RuntimeError('freq_downsamp does not evenly divide NCHAN ({0}).'.format(nchan))

        nrows_total = end_row - start_row + 1
        down_nsamp_blk = nsamp_blk / downsamp
        down_nchan = nchan / freq_downsamp
        down_tbin = tbin * downsamp
        # print('nrows_total: {0}'.format(nrows_total))
        # print('down_nsamp_blk: {0}'.format(down_nsamp_blk))
        # print('down_nchan: {0}'.format(down_nchan))
        # print('down_tbin: {0}'.format(down_tbin))

        # Data types of the signed and unsigned.
        if nbit == 8:
            signed_type = numpy.int8
            unsigned_type = numpy.uint8
        elif nbit == 16:
            signed_type = numpy.int16
            unsigned_type = numpy.uint16
        elif nbit == 32:
            signed_type = numpy.float32
            unsigned_type = numpy.float32
        else:
            raise RuntimeError('Not handled number of NBITS ({0}).'.format(nbit))

        # Allocate arrays.
        row_spectrum_result = numpy.zeros(nchan, dtype=numpy.float32)
        final_spectrum = numpy.zeros((nrows_total * down_nsamp_blk, npol, down_nchan), dtype=numpy.float32)
        if get_ft:
            final_freqs = numpy.zeros(down_nchan)
            final_timing = numpy.zeros(nrows_total * down_nsamp_blk)
        pol_sign = 1
        if 'AABB' in pol_type:
            pol_sign = 2

        # Iterate over rows (blocks).
        for irow in range(nrows_total):
            if apply_scales:
                offsets = self.fits['SUBINT'].data['DAT_OFFS'][irow + start_row]
                scales = self.fits['SUBINT'].data['DAT_SCL'][irow + start_row]
                weights = self.fits['SUBINT'].data['DAT_WTS'][irow + start_row]
                offsets = offsets.reshape((npol, nchan))
                scales = scales.reshape((npol, nchan))
                # weights = numpy.concatenate((weights, weights, weights, weights))
                # weights = weights.reshape((npol, nchan))
            if get_ft:
                row_timing = self.fits['SUBINT'].data['OFFS_SUB'][irow + start_row] - (self.fits['SUBINT'].data['TSUBINT'][irow + start_row] / 2.0)
                freqs_row = self.fits['SUBINT'].data['DAT_FREQ'][irow + start_row]
            row_data = self.fits['SUBINT'].data['DATA'][irow + start_row]
            # print('row_data: {0}'.format(row_data.shape))
            # print('offets: {0}'.format(offsets.shape))
            # print('scales: {0}'.format(scales.shape))
            # print('weights: {0}'.format(weights.shape))

        # Decode row_data for 16 bit data type file.
        if (nbit == 16):
            row_data = numpy.fromstring(row_data.tostring(), dtype=numpy.int16)
            row_data = row_data.reshape((nsamp_blk, npol, nchan))

        # Iterate over samples and polarisations.
        for isamp in range(down_nsamp_blk):
            for ipol in range(npol):
                if ipol < pol_sign:
                    data_type = unsigned_type
                else:
                    data_type = signed_type
                row_spectrum_result = row_data[isamp * downsamp:(isamp + 1) * downsamp, ipol, :].astype(data_type).mean(0)
                if get_ft:
                    final_timing[irow * down_nsamp_blk + isamp] = row_timing + (isamp + 0.5) * down_tbin
                if apply_scales:
                    row_spectrum_result *= scales[ipol, :]
                    row_spectrum_result += offsets[ipol, :]
                    row_spectrum_result *= weights[ipol, :]
                if freq_downsamp == 1:
                    final_spectrum[irow * down_nsamp_blk + isamp, ipol, :] = row_spectrum_result
                    if get_ft:
                        final_freqs[:] = freqs_row[:]
                else:
                    final_spectrum[irow * down_nsamp_blk + isamp, ipol, :] = row_spectrum_result.reshape((-1, freq_downsamp)).mean(1)
                    if get_ft:
                        final_freqs[:] = freqs_row.reshape((-1, freq_downsamp)).mean(1)

        if squeeze:
            final_spectrum = final_spectrum.squeeze()
        if get_ft:
            return (final_spectrum, final_timing, final_freqs)
        else:
            return final_spectrum
