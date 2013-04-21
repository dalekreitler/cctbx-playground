from __future__ import division

# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id: cspad_detector.py 371 2013-04-14 18:25:42Z sauter $

import math

from iotbx.detectors import generic_flex_image
from libtbx import easy_pickle
from scitbx.array_family import flex
from scitbx.matrix import col, rec, sqr
from xfel.cftbx.detector.generic_detector import GenericDetector

class CSPadDetector(GenericDetector):
  def __init__(self, filename):
    self.filename = filename
    self.pixel_size = 0.1
    self.vendortype = "npy_raw"
    self.beamx = 0.
    self.beamy = 0.

  def __getattr__(self, attr):
    if   attr=='attenuation' : return self._metrology_params.attenuation
    elif attr=='beam_center' : return self._metrology_params.beam_center
    elif attr=='distance' : return self._metrology_params.distance
    elif attr=='pulse_length' : return self._metrology_params.pulse_length
    elif attr=='sequence_number' : return self._metrology_params.sequence_number
    elif attr=='timestamp' : return self._metrology_params.timestamp
    elif attr=='wavelength' : return self._metrology_params.wavelength
    elif attr=='xtal_target' : return self._metrology_params.xtal_target
    else: raise AttributeError(attr)

  def readHeader(self):
    from xfel.cftbx.detector.metrology import metrology_as_transformation_matrices

    d = easy_pickle.load(self.filename)

    # Derive the transformation matrices from the metrology in the
    # image.
    self._metrology_params = d["METROLOGY"].extract()
    self._matrices = metrology_as_transformation_matrices(
      self._metrology_params)

    self._tiles = d["TILES"]
    self._keylist = self._tiles.keys() #for later use by get_pixel_intensity()

    # Assert that all ASIC:s are the same size, and that there are
    # transformation matrices for each ASIC.
    for (key, asic) in self._tiles.iteritems():
      if not hasattr(self, "_asic_focus"):
        self._asic_focus = asic.focus()
      else:
        assert asic.focus() == self._asic_focus
      assert key in self._matrices

    # Assert that all specified pixel sizes and saturated values are
    # equal.
    for p in self._metrology_params.detector.panel:
      for s in p.sensor:
        for a in s.asic:
          if not hasattr(self, "_pixel_size"):
            self._pixel_size = a.pixel_size
          else:
            assert self._pixel_size == a.pixel_size

          if not hasattr(self, "_saturation"):
            self._saturation = a.saturation
          else:
            # XXX real-valued equality!  See
            # cctbx_project/scitbx/math/approx_equal.h
            assert self._saturation == a.saturation
    assert hasattr(self, "_pixel_size")
    assert hasattr(self, "_saturation")

    # Determine next multiple of eight.  Set size1 and size2 to the
    # focus of the padded rawdata.
    self._asic_padded = (8 * int(math.ceil(self._asic_focus[0] / 8)),
                         8 * int(math.ceil(self._asic_focus[1] / 8)))
    self.size1 = len(self._tiles) * self._asic_padded[0]
    self.size2 = self._asic_padded[1]


  def show_header(self):
    return "CSPad detector with nothing in it"


  def read(self):
    pass

  def supports_quadrant_calibration(self): return True
  def implements_metrology_matrices(self): return True

  def apply_metrology_from_matrices(self, matrices):
    """The apply_metrology_from_matrices() function replaces the
    current set of transformation matrices with that given in the
    dictionary @p matrices.

    XXX Could implement metrology "adjustment" as opposed metrology
    "replacement".
    """

    self._matrices = matrices

    # Old apply_metrology_from_params() implementation, bound to be
    # useful somewhere.
    """
    from xfel.cftbx.detector.metrology import metrology_as_transformation_matrices
    from libtbx.phil import experimental

    experimental.merge_params_by_key(self._metrology_params, params, "serial")
    self._matrices = metrology_as_transformation_matrices(
      self._metrology_params)

#    from xfel.cftbx.detector.metrology2phil import master_phil
#    master_phil.format(python_object=self._metrology_params).show()
    """


  def get_flex_image(self, brightness, **kwargs):
    # no kwargs supported at present
    """XXX"""

    # Horizontal and vertical dimensions of a pixel, in meters.  XXX
    # Sort out relationship to pixel_size member!
    pixel_size = [110e-6, 110e-6]

    # E maps picture coordinates onto metric Cartesian coordinates,
    # i.e. [row, column, 1 ] -> [x, y, z, 1].  Both frames share the
    # same origin, but the first coordinate of the screen coordinate
    # system increases downwards, while the second increases towards
    # the right.  XXX Is this orthographic projection the only one
    # that makes any sense?
    E = rec(elems=[0, +pixel_size[1], 0,
                   -pixel_size[0], 0, 0,
                   0, 0, 0,
                   0, 0, 1],
            n=[4, 3])

    # P: [x, y, z, 1] -> [row, column, 1] XXX Projection matrix should
    # probably be implemented elsewhere.
    P = rec(
      elems=[0, -1 / pixel_size[1], 0, (self._asic_focus[0] - 1) / 2,
             +1 / pixel_size[0], 0, 0, (self._asic_focus[1] - 1) / 2,
             0, 0, 0, 1],
      n=[3, 4])

    # XXX Add ASIC:s in order?  If a point is contained in two ASIC:s
    # simultaneously, it will be assigned to the ASIC defined first.
    # XXX Use a Z-buffer instead?
    nmemb = 0
    for key, asic in self._tiles.iteritems():
      # Create my_flex_image and rawdata on the first iteration.
      if ("rawdata" not in locals()):
        rawdata = flex.double(flex.grid(self.size1, self.size2))
        my_flex_image = generic_flex_image(
          rawdata=rawdata,
          size1_readout=self._asic_focus[0],
          size2_readout=self._asic_focus[1],
          brightness=brightness,
          saturation=self._saturation)

      rawdata.matrix_paste_block_in_place(
        block=asic,
        i_row=nmemb * self._asic_padded[0],
        i_column=0)
      nmemb += 1

      # key is guaranteed to exist in self._matrices as per
      # readHeader().  Last row of self._matrices[key][0] is always
      # [0, 0, 0, 1].
      T = P * self._matrices[key][0] * E
      R = sqr([T(0, 0), T(0, 1),
               T(1, 0), T(1, 1)])
      t = col([T(0, 2), T(1, 2)])

      my_flex_image.add_transformation_and_translation(R, t)
    my_flex_image.followup_brightness_scale()
    return my_flex_image


  def get_panel_fast_slow(self, serial):
    """Get the average x- and y-coordinates of all the ASIC:s in the
    panel with serial @p serial.  This is done by back-transforming
    the centre's of each ASIC to the screen (sort of) coordinate
    system.  This is more robust than getting the panel positions
    directly.
    """

    center = col([self._asic_focus[0] / 2, self._asic_focus[1] / 2, 1])
    fast, nmemb, slow = 0, 0, 0

    # Use the pixel size for the ASIC to construct the final
    # projection matrix.
    for p in self._metrology_params.detector.panel:
      if (p.serial != serial):
        continue
      for s in p.sensor:
        for a in s.asic:
          E = rec(elems=[+1 / a.pixel_size[0], 0, 0, 0,
                         0, -1 / a.pixel_size[1], 0, 0],
                  n=[2, 4])
          P = rec(elems=[
              0, +a.pixel_size[0], a.pixel_size[0] * (1 - a.dimension[0]) / 2,
              -a.pixel_size[1], 0, a.pixel_size[1] * (a.dimension[1] - 1) / 2,
              0, 0, 0,
              0, 0, 1],
                  n=[4, 3])
          T = self._matrices[(0, p.serial, s.serial, a.serial)][1]

          t = E * T * P * center
          fast += t(0, 0)
          slow += t(1, 0)
          nmemb += 1
    if (nmemb == 0):
      return (0, 0)
    return (fast / nmemb, slow / nmemb)


  def displace_panel_fast_slow(self, serial, fast, slow):
    """Displace all ASICS:s in the panel with serial @p serial such
    that their new average position becomes @p fast, @p slow.  The
    function returns the updated transformation matrices.
    """

    # XXX Get from elsewhere!
    pixel_size = [110e-6, 110e-6]

    dx, dy = fast * pixel_size[0], -slow * pixel_size[1]

    for key, (Tf, Tb) in self._matrices.iteritems():
      if (len(key) == 4 and key[1] == serial):
        Tb_new = sqr(
          [Tb(0, 0), Tb(0, 1), Tb(0, 2), Tb(0, 3) + dx,
           Tb(1, 0), Tb(1, 1), Tb(1, 2), Tb(1, 3) + dy,
           Tb(2, 0), Tb(2, 1), Tb(2, 2), Tb(2, 3) + 0,
           Tb(3, 0), Tb(3, 1), Tb(3, 2), Tb(3, 3) + 0])

        # XXX Math worked out elsewhere.
        Tf_new = sqr(
          [Tf(0, 0), Tf(0, 1), Tf(0, 2), Tf(0, 3) - Tf(0, 0) * dx - Tf(0, 1) * dy,
           Tf(1, 0), Tf(1, 1), Tf(1, 2), Tf(1, 3) - Tf(1, 0) * dx - Tf(1, 1) * dy,
           Tf(2, 0), Tf(2, 1), Tf(2, 2), Tf(2, 3) - Tf(2, 0) * dx - Tf(2, 1) * dy,
           Tf(3, 0), Tf(3, 1), Tf(3, 2), Tf(3, 3) - Tf(3, 0) * dx - Tf(3, 1) * dy])

        self._matrices[key] = (Tf_new, Tb_new)
    return self._matrices


  def get_pixel_intensity(self,coords):
    tileno = int(coords[2])
    if tileno < 0: return None
    try:
      return self._tiles[self._keylist[tileno]][
        (int(round(coords[0],0)), (int(round(coords[1],0))))]
    except IndexError:
      return None


  def _average_transformation(self, matrices, keys):
    """The _average_transformation() function determines the average
    rotation and translation from the transformation matrices in @p
    matrices with keys matching @p keys.  The function returns a
    two-tuple of the average rotation in quaternion representation and
    the average translation.
    """

    from xfel.cftbx.detector.metrology import _transform # XXX private!
    from scitbx.matrix import zeros
    from tntbx import svd

    # Sum all rotation matrices and translation vectors.
    sum_R = flex.double(flex.grid(3, 3))
    sum_t = flex.double(flex.grid(3, 1))
    nmemb = 0
    for key in keys:
      T = matrices[key][1]
      sum_R += flex.double([T(0, 0), T(0, 1), T(0, 2),
                            T(1, 0), T(1, 1), T(1, 2),
                            T(2, 0), T(2, 1), T(2, 2)])
      sum_t += flex.double([T(0, 3), T(1, 3), T(2, 3)])
      nmemb += 1
    if (nmemb == 0):
      # Return zero-rotation and zero-translation.
      return (col([1, 0, 0, 0]), zeros((3, 1)))

    # Calulate average rotation matrix as U * V^T where sum_R = U * S
    # * V^T and S diagonal (Curtis et al. (1993) 377-385 XXX proper
    # citation, repeat search), and convert to quaternion.
    svd = svd(sum_R)
    R_avg = sqr(list(svd.u().matrix_multiply_transpose(svd.v())))
    o_avg = R_avg.r3_rotation_matrix_as_unit_quaternion()
    t_avg = col(list(sum_t / nmemb))

    return (o_avg, t_avg)


  def _regularize_transformations(self, matrices, key=(0,)):
    """The _regularize_transformations() function recursively sets the
    transformation of each element to the average transformation of
    its children.  XXX Round average orientations to multiples of 90
    degrees?  XXX Should they maybe return matrices instead?.
    """

    from xfel.cftbx.detector.metrology import _transform # XXX private!

    # The length of the key tuple identifies the current recursion
    # depth.
    next_depth = len(key) + 1

    # Base case: return if at the maximum depth.
    keys = [k for k in matrices.keys()
            if (len(k) == next_depth and k[0:next_depth - 1] == key)]
    if (len(keys) == 0):
      return

    # Recursion: regularize at next depth, and get the average.
    for k in keys:
      self._regularize_transformations(matrices, k)
    o, t = self._average_transformation(matrices, keys)
    matrices[key] = _transform(o, t)


  def _matrix_as_string(self, T):
    """XXX other uses?"""

    R = sqr([T(0, 0), T(0, 1), T(0, 2),
             T(1, 0), T(1, 1), T(1, 2),
             T(2, 0), T(2, 1), T(2, 2)])
    o = R.r3_rotation_matrix_as_unit_quaternion()
    t = col([T(0, 3), T(1, 3), T(2, 3)])

    return ("orientation = %s, %s, %s, %s\n" % tuple(repr(c) for c in o)) + \
        ("translation = %s, %s, %s\n" % tuple(repr(c) for c in t))


  def transformation_matrices_as_metrology(self):
    """The transformation_matrices_as_metrology() function regularizes
    the the transformation matrices and converts them to a phil
    object.
    """

    from xfel.cftbx.detector.metrology2phil import master_phil
    from libtbx import phil

    # XXX Experimental!
    self._regularize_transformations(self._matrices)

    (Tf_d, Tb_d) = self._matrices[(0,)]
    metrology_str = "detector {\n"
    metrology_str += "serial = %d\n" % 0
    metrology_str += self._matrix_as_string(Tb_d)

    for p in [k[1] for k in self._matrices.keys()
              if (len(k) == 2 and k[0:1] == (0,))]:
      (Tf_p, Tb_p) = self._matrices[(0, p)]
      metrology_str += "panel {\n"
      metrology_str += "serial = %d\n" % p
      metrology_str += self._matrix_as_string(Tf_d * Tb_p)

      for s in[k[2] for k in self._matrices.keys()
               if (len(k) == 3 and k[0:2] == (0, p))]:
        (Tf_s, Tb_s) = self._matrices[(0, p, s)]
        metrology_str += "sensor {\n"
        metrology_str += "serial = %d\n" % s
        metrology_str += self._matrix_as_string(Tf_p * Tb_s)

        for a in [k[3] for k in self._matrices.keys()
                  if (len(k) == 4 and k[0:3] == (0, p, s))]:
          (Tf_a, Tb_a) = self._matrices[(0, p, s, a)]
          metrology_str += "asic {\n"
          metrology_str += "serial = %d\n" % a
          metrology_str += self._matrix_as_string(Tf_s * Tb_a)

          metrology_str += "}\n"
        metrology_str += "}\n"
      metrology_str += "}\n"
    metrology_str += "}\n"

    return master_phil.fetch(sources=[phil.parse(metrology_str)])


  def readout_coords_as_detector_coords(self, coords):
    """
    Convert a 3 tuple coordinates from readout space (x, y, tile number)
    to detector space in meters, relative to the detector center
    """
    tileno = int(coords[2])
    if tileno < 0: return None
    try:
      T_f, T_b = self._matrices[self._keylist[tileno]]
    except IndexError:
      return None
    assert self._pixel_size is not None
    assert self._asic_focus is not None

    from xfel.cftbx.detector.metrology import get_projection_matrix

    P_f, P_b = get_projection_matrix(self._pixel_size, self._asic_focus)

    return T_b * P_b * col([int(coords[0]), int(coords[1]), 1])

  def image_coords_as_detector_coords (self, x, y, readout=None) :
      """
      Convert image pixel coordinates to absolute position on the detector
      (in mm).
      """
      if readout is None:
        return super(CSPadDetector,self).image_coords_as_detector_coords(x, y)

      c = self.readout_coords_as_detector_coords([x, y, readout])
      return c[0] * 1000, c[1] * 1000

  def bounding_box_mm (self) :
      """
      Calculate the extent of this tiled detector image in mm.
      Returns x , y, width, height, where x and y are the upper left coordinates of the detector
      (which may or may not actually have a pixel if the nearest asic is tilted)
      """
      left = top = float("inf")
      right = bottom = -float("inf")
      for tileno in range(0,len(self._tiles)):
          coords = self.readout_coords_as_detector_coords((0,0,tileno))
          if(coords[0] < left): left = coords[0]
          elif(coords[0] > right): right = coords[0]
          if(coords[1] < top): top = coords[1]
          elif(coords[1] > bottom): bottom = coords[1]

          coords = self.readout_coords_as_detector_coords((self._asic_focus[0]-1,0,tileno))
          if(coords[0] < left): left = coords[0]
          elif(coords[0] > right): right = coords[0]
          if(coords[1] < top): top = coords[1]
          elif(coords[1] > bottom): bottom = coords[1]

          coords = self.readout_coords_as_detector_coords((0,self._asic_focus[1]-1,tileno))
          if(coords[0] < left): left = coords[0]
          elif(coords[0] > right): right = coords[0]
          if(coords[1] < top): top = coords[1]
          elif(coords[1] > bottom): bottom = coords[1]

          coords = self.readout_coords_as_detector_coords((self._asic_focus[0]-1,self._asic_focus[1]-1,tileno))
          if(coords[0] < left): left = coords[0]
          elif(coords[0] > right): right = coords[0]
          if(coords[1] < top): top = coords[1]
          elif(coords[1] > bottom): bottom = coords[1]

      return (left*1000, top*1000, (right-left)*1000, (bottom-top)*1000)

  def detector_coords_as_image_coords_float (self, x, y) :
    """
    Convert absolute detector position (in mm) to floating-value image pixel coordinates.
    """
    if self._pixel_size is None or type(self._pixel_size) is float:
      return super(CSPadDetector,self).detector_coords_as_image_coords_float(x, y)

    return x / self._pixel_size[0] / 1000, \
           y / self._pixel_size[1] / 1000

  def get_raw_data(self):
    # Not intended for production; simply a means to marshall all same-size tile
    # data together to report it out as a single array; used for testing dxtbx.
    keys = self._tiles.keys()
    keys.sort()
    raw = flex.double(flex.grid(len(keys)*self._tiles[keys[0]].focus()[0],
                                          self._tiles[keys[0]].focus()[1]))
    slowstride = self._tiles[keys[0]].focus()[0]
    for ik,k in enumerate(keys):
      raw.matrix_paste_block_in_place(
        block = self._tiles[k],
        i_row = slowstride * ik,
        i_column = 0)
    return raw
