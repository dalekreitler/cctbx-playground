from __future__ import division
#!/usr/bin/env python
# scan.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# A model for the scan for the "updated experimental model" project documented
# in internal ticket #1555. This is not designed to be used outside of the
# XSweep classes.

import pycbf
import copy
from dxtbx_model_ext import Scan

from scan_helpers import scan_helper_image_files
from scan_helpers import scan_helper_image_formats


class scan_factory:
  '''A factory for scan instances, to help with constructing the classes
  in a set of common circumstances.'''

  @staticmethod
  def make_scan(image_range, exposure_times, oscillation, epochs, deg=True):
    from scitbx.array_family import flex

    if not isinstance(exposure_times, list):
      num_images = image_range[1] - image_range[0] + 1
      exposure_times = [exposure_times for i in range(num_images)]
    else:
      num_images = image_range[1] - image_range[0] + 1
      num_exp = len(exposure_times)
      if num_exp != num_images:
        if (num_exp == 0):
          exposure_times = [0 for i in range(num_images)]
        else:
          exposure_times = exposure_times.extend(
            [exposure_times[-1] for i in range(num_images - num_exp)])

    epoch_list = [epochs[j] for j in sorted(epochs)]

    return Scan(
        tuple(map(int, image_range)),
        tuple(map(float, oscillation)),
        flex.double(list(map(float, exposure_times))),
        flex.double(list(map(float, epoch_list))),
        deg)

  @staticmethod
  def single(filename, format, exposure_times, osc_start, osc_width, epoch):
    '''Construct an scan instance for a single image.'''

    import os
    index = scan_helper_image_files.image_to_index(os.path.split(filename)[-1])
    if epoch is None:
      epoch = 0.0
    return scan_factory.make_scan(
                (index, index), exposure_times, (osc_start, osc_width),
                {index:epoch})

  @staticmethod
  def imgCIF(cif_file):
    '''Initialize a scan model from an imgCIF file.'''

    cbf_handle = pycbf.cbf_handle_struct()
    cbf_handle.read_file(cif_file, pycbf.MSG_DIGEST)

    return scan_factory.imgCIF_H(cif_file, cbf_handle)

  @staticmethod
  def imgCIF_H(cif_file, cbf_handle):
    '''Initialize a scan model from an imgCIF file handle, where it is
    assumed that the file has already been read.'''

    exposure = cbf_handle.get_integration_time()
    timestamp = cbf_handle.get_timestamp()[0]

    gonio = cbf_handle.construct_goniometer()
    angles = tuple(gonio.get_rotation_range())

    # xia2-56 handle gracefully reverse turning goniometers - this assumes the
    # rotation axis is correctly inverted in the goniometer factory
    if angles[1] < 0:
      angles = -angles[0], -angles[1]

    index = scan_helper_image_files.image_to_index(cif_file)

    gonio.__swig_destroy__(gonio)

    return scan_factory.make_scan(
        (index, index), exposure, angles, {index:timestamp})

  @staticmethod
  def add(scans):
    '''Sum a list of scans wrapping the sligtly clumsy idiomatic method:
    sum(scans[1:], scans[0]).'''
    return sum(scans[1:], scans[0])

  @staticmethod
  def search(filename):
    '''Get a list of files which appear to match the template and
    directory implied by the input filename. This could well be used
    to get a list of image headers to read and hence construct scans
    from.'''

    template, directory = \
              scan_helper_image_files.image_to_template_directory(filename)

    indices = scan_helper_image_files.template_directory_to_indices(
        template, directory)

    return [scan_helper_image_files.template_directory_index_to_image(
        template, directory, index) for index in indices]

  @staticmethod
  def format(name):
    '''Return the correct format token for a given name, for example:

    cbf, CBF
    smv, SMV
    tiff, tif, TIFF
    raxis, RAXIS
    mar, MAR

    to the appropriate static token which will be used as a handle
    everywhere else in this.'''

    if name.upper() == 'CBF':
      return scan_helper_image_formats.FORMAT_CBF
    elif name.upper() == 'SMV':
      return scan_helper_image_formats.FORMAT_SMV
    elif name.upper() == 'TIF' or name.upper() == 'TIFF':
      return scan_helper_image_formats.FORMAT_TIFF
    elif name.upper() == 'RAXIS':
      return scan_helper_image_formats.FORMAT_RAXIS
    elif name.upper() == 'MAR':
      return scan_helper_image_formats.FORMAT_MAR

    raise RuntimeError, 'name %s not known' % name
