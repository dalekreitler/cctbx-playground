from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.map_comparison

from cctbx import crystal
from cctbx import maptbx, miller
from cctbx.sgtbx import space_group_info
from iotbx import file_reader, phil
import iotbx.ccp4_map
from libtbx.utils import Sorry
from scitbx.array_family import flex
import os, sys

master_phil = phil.parse("""
include scope libtbx.phil.interface.tracking_params
input
{
  map_1 = None
    .type = path
    .short_caption = Map 1
    .help = A CCP4-formatted map
    .style = file_type:ccp4_map bold input_file
  map_2 = None
    .type = path
    .short_caption = Map 2
    .help = A CCP4-formatted map
    .style = file_type:ccp4_map bold input_file
  mtz_1 = None
    .type = path
    .short_caption = Map 1
    .help = MTZ file containing map
    .style = file_type:hkl bold input_file
  mtz_2 = None
    .type = path
    .short_caption = Map 2
    .help = MTZ file containing map
    .style = file_type:hkl bold input_file
  mtz_label_1 = None
    .type = str
    .short_caption = Data label
    .help = Data label for complex map coefficients in MTZ file
  mtz_label_2 = None
    .type = str
    .short_caption = Data label
    .help = Data label for complex map coefficients in MTZ file
}
options
{
  resolution_factor = 0.25
    .type = float
    .short_caption = Resolution factor
    .help = Determines grid spacing in map
}
""", process_includes=True)

master_params = master_phil

def show_overall_statistics(out=sys.stdout, s=None, header=None):
  print >> out, header
  print >> out, "  min/max/mean: %6.4f %6.4f %6.4f"%(s.min(), s.max(), s.mean())
  print >> out, "  kurtosis    : %6.4f" % s.kurtosis()
  print >> out, "  skewness    : %6.4f" % s.skewness()
  print >> out, "  sigma       : %6.4f" % s.sigma()

def create_statistics_dict(out=sys.stdout, s=None):
  statistics_dict = dict()
  statistics_dict['min'] = s.min()
  statistics_dict['max'] = s.max()
  statistics_dict['mean'] = s.mean()
  statistics_dict['kurtosis'] = s.kurtosis()
  statistics_dict['skewness'] = s.skewness()
  statistics_dict['sigma'] = s.sigma()
  return statistics_dict

def show_citation(out=sys.stdout):
  print >> out, "-"*79
  msg = """Map comparison and statistics. For details see:
  Acta Cryst. (2014). D70, 2593-2606
  Metrics for comparison of crystallographic maps
  A. Urzhumtsev, P. V. Afonine, V. Y. Lunin, T. C. Terwilliger and P. D. Adams"""
  print >> out, msg
  print >> out, "-"*79

# =============================================================================
def run(args, out=sys.stdout, validated=False):
  show_citation(out=out)
  if (len(args) == 0):
    master_phil.show(out=out)
    print >> out,\
      '\nUsage: phenix.map_comparison <CCP4> <CCP4>\n',\
      '       phenix.map_comparison <CCP4> <MTZ> mtz_label_1=<label>\n',\
      '       phenix.map_comparison <MTZ 1> mtz_label_1=<label 1> <MTZ 2> mtz_label_2=<label 2>\n'
    sys.exit()

  # process arguments
  params = None
  input_attributes = ['map_1', 'map_2', 'mtz_1', 'mtz_2']
  try: # automatic parsing
    params = phil.process_command_line_with_files(
      args=args, master_phil=master_phil).work.extract()
  except Exception: # map_file_def only handles one map phil
    from libtbx.phil.command_line import argument_interpreter
    arg_int = argument_interpreter(master_phil=master_phil)
    command_line_args = list()
    map_files = list()
    for arg in args:
      if (os.path.isfile(arg)):
        map_files.append(arg)
      else:
        command_line_args.append(arg_int.process(arg))
    params = master_phil.fetch(sources=command_line_args).extract()

    # check if more files are necessary
    n_defined = 0
    for attribute in input_attributes:
      if (getattr(params.input, attribute) is not None):
        n_defined += 1

    # matches files to phil scope, stops once there is sufficient data
    for map_file in map_files:
      if (n_defined < 2):
        current_map = file_reader.any_file(map_file)
        if (current_map.file_type == 'ccp4_map'):
          n_defined += 1
          if (params.input.map_1 is None):
            params.input.map_1 = map_file
          elif (params.input.map_2 is None):
            params.input.map_2 = map_file
        elif (current_map.file_type == 'hkl'):
          n_defined += 1
          if (params.input.mtz_1 is None):
            params.input.mtz_1 = map_file
          elif (params.input.mtz_2 is None):
            params.input.mtz_2 = map_file
      else:
        print >> out, 'WARNING: only the first two files are used'
        break

  # validate arguments (GUI sets validated to true, no need to run again)
  assert (params is not None)
  if (not validated):
    validate_params(params)

  # ---------------------------------------------------------------------------
  # check if maps need to be generated from mtz
  n_maps = 0
  maps = list()
  for attribute in input_attributes:
    filename = getattr(params.input, attribute)
    if (filename is not None):
      current_map = file_reader.any_file(filename)
      maps.append(current_map)
      if (current_map.file_type == 'ccp4_map'):
        n_maps += 1

  # construct maps, if necessary
  crystal_gridding = None
  m1 = None
  m2 = None

  # 1 map, 1 mtz file
  if (n_maps == 1):
    for current_map in maps:
      if (current_map.file_type == 'ccp4_map'):
        uc = current_map.file_object.unit_cell()
        sg_info = space_group_info(current_map.file_object.space_group_number)
        n_real = current_map.file_object.unit_cell_grid
        crystal_gridding = maptbx.crystal_gridding(
          uc, space_group_info=sg_info, pre_determined_n_real=n_real)
        m1 = current_map.file_object.map_data()
    if (crystal_gridding is not None):
      label = None
      for attribute in [('mtz_1', 'mtz_label_1'),
                        ('mtz_2', 'mtz_label_2')]:
        filename = getattr(params.input, attribute[0])
        label = getattr(params.input, attribute[1])
        if ( (filename is not None) and (label is not None) ):
          break
      # labels will match currently open mtz file
      for current_map in maps:
        if (current_map.file_type == 'hkl'):
          m2 = miller.fft_map(
            crystal_gridding=crystal_gridding,
            fourier_coefficients=current_map.file_server.get_miller_array(
              label)).apply_sigma_scaling().real_map_unpadded()
    else:
      raise Sorry('Gridding is not defined.')

  # 2 mtz files
  elif (n_maps == 0):
    crystal_symmetry = get_crystal_symmetry(maps[0])
    d_min = min(get_d_min(maps[0]), get_d_min(maps[1]))
    crystal_gridding = maptbx.crystal_gridding(
      crystal_symmetry.unit_cell(), d_min=d_min,
      resolution_factor=params.options.resolution_factor,
      space_group_info=crystal_symmetry.space_group_info())
    m1 = miller.fft_map(
      crystal_gridding=crystal_gridding,
      fourier_coefficients=maps[0].file_server.get_miller_array(
        params.input.mtz_label_1)).apply_sigma_scaling.real_map_unpadded()
    m2 = miller.fft_map(
      crystal_gridding=crystal_gridding,
      fourier_coefficients=maps[1].file_server.get_miller_array(
        params.input.mtz_label_2)).apply_sigma_scaling.real_map_unpadded()

  # 2 maps
  else:
    m1 = maps[0].file_object.map_data()
    m2 = maps[1].file_object.map_data()

  # ---------------------------------------------------------------------------
  # analyze maps
  assert ( (m1 is not None) and (m2 is not None) )

  # show general statistics
  s1 = maptbx.more_statistics(m1)
  s2 = maptbx.more_statistics(m2)
  show_overall_statistics(s=s1, header="Map 1 (%s):"%params.input.map_1)
  show_overall_statistics(s=s2, header="Map 2 (%s):"%params.input.map_2)
  cc_input_maps = flex.linear_correlation(x = m1.as_1d(),
                                          y = m2.as_1d()).coefficient()
  print >> out, "CC, input maps: %6.4f" % cc_input_maps

  # compute CCpeak
  cc_peaks = list()
  m1_he = maptbx.volume_scale(map = m1,  n_bins = 10000).map_data()
  m2_he = maptbx.volume_scale(map = m2,  n_bins = 10000).map_data()
  cc_quantile = flex.linear_correlation(x = m1_he.as_1d(),
                                        y = m2_he.as_1d()).coefficient()
  print >> out, "CC, quantile rank-scaled (histogram equalized) maps: %6.4f" % \
    cc_quantile
  print >> out, "Peak correlation:"
  print >> out, "  cutoff  CCpeak"
  for cutoff in [i/100. for i in range(0,100,5)]+[0.99, 1.0]:
    cc_peak = maptbx.cc_peak(map_1=m1_he, map_2=m2_he, cutoff=cutoff)
    print >> out, "  %3.2f   %7.4f" % (cutoff, cc_peak)
    cc_peaks.append((cutoff, cc_peak))

  # compute discrepancy function (D-function)
  discrepancies = list()
  cutoffs = flex.double([i/20. for i in range(1,20)])
  df = maptbx.discrepancy_function(map_1=m1_he, map_2=m2_he, cutoffs=cutoffs)
  print >> out, "Discrepancy function:"
  print >> out, "  cutoff  D"
  for c, d in zip(cutoffs, df):
    print >> out, "  %3.2f   %7.4f" % (c,d)
    discrepancies.append((c, d))

  # compute and output histograms
  h1 = maptbx.histogram(map=m1, n_bins=10000)
  h2 = maptbx.histogram(map=m2, n_bins=10000)
  print >> out, "Map histograms:"
  print >> out, "Map 1 (%s)     Map 2 (%s)"%\
    (params.input.map_1,params.input.map_2)
  print >> out, "(map_value,cdf,frequency) <> (map_value,cdf,frequency)"
  for a1,c1,v1, a2,c2,v2 in zip(h1.arguments(), h1.c_values(), h1.values(),
                                h2.arguments(), h2.c_values(), h2.values()):
    print >> out, "(%9.5f %9.5f %9.5f) <> (%9.5f %9.5f %9.5f)"%\
      (a1,c1,v1, a2,c2,v2)

  # store results
  s1_dict = create_statistics_dict(s=s1)
  s2_dict = create_statistics_dict(s=s2)
  results = dict()
  results['map_files'] = (params.input.map_1, params.input.map_2)
  results['map_statistics'] = (s1_dict, s2_dict)
  results['cc_input_maps'] = cc_input_maps
  results['cc_quantile'] = cc_quantile
  results['cc_peaks'] = cc_peaks
  results['discrepancies'] = discrepancies
  results['map_histograms'] = ( (h1.arguments(), h1.c_values(), h1.values()),
                                (h2.arguments(), h2.c_values(), h2.values()) )

  return results

# -----------------------------------------------------------------------------
def get_crystal_symmetry(file_handle):
  '''
  Helper function for get crystal symmetry from files
  '''
  file_object = file_handle.file_object
  cs = None
  if (hasattr(file_object, 'space_group_number')):     # CCP4 map
    cs = crystal.symmetry(file_object.unit_cell().parameters(),
                          file_object.space_group_number)
  elif (hasattr(file_object, 'as_miller_arrays')):     # MTZ file
    ma = file_object.as_miller_arrays()
    for a in ma:
      if (a.is_complex_array()):
        cs = a.crystal_symmetry()
        break
  if (cs is None):
    raise Sorry('Could not find crystal symmetry in %s.' %
                file_handle.file_name)
  return cs

def get_d_min(file_handle):
  '''
  Helper function for getting d_min from mtz file
  '''
  miller_arrays = file_handle.file_server.miller_arrays
  d_min = 10.0
  for miller_array in miller_arrays:
    d_min = min(d_min, miller_array.d_min())
  return d_min

def get_mtz_labels(file_handle):
  '''
  Helper function for getting data labels for complex miller arrays
  Returns list of labels for complex arrays
  '''
  miller_arrays = file_handle.file_server.miller_arrays
  labels = list()
  for miller_array in miller_arrays:
    if (miller_array.is_complex_array()):
      labels.append(miller_array.info().label_string())
  return labels

# =============================================================================
# Parameter validation for CLI and GUI
def validate_params(params):

  # check that only 2 files, in any combination, are provided
  input_attributes = ['map_1', 'map_2', 'mtz_1', 'mtz_2']
  n_defined = 0
  for attribute in input_attributes:
    if (getattr(params.input, attribute) is not None):
      n_defined += 1
  if (n_defined != 2):
    raise Sorry('Insufficient data, please provide 2 files' +
                ' (CCP4-formated map or MTZ)')

  # check file type
  maps = list()
  for attribute in input_attributes:
    filename = getattr(params.input, attribute)
    if (filename is not None):
      file_handle = file_reader.any_file(filename)
      if ( (file_handle.file_type != 'ccp4_map') and
           (file_handle.file_type != 'hkl') ):
        raise Sorry('Please input a CCP4-formatted map or MTZ file for %s.'\
                    % filename)
      else:
        maps.append(file_handle)

  # check symmetry
  cs1 = get_crystal_symmetry(maps[0])
  cs2 = get_crystal_symmetry(maps[1])
  if (cs1.is_similar_symmetry(cs2) is False):
    raise Sorry('The symmetry of the two files is not similar.')

  # check gridding if 2 map files are provided
  if ( (maps[0].file_type == 'ccp4_map') and
       (maps[1].file_type == 'ccp4_map') ):
    m1 = maps[0].file_object.map_data()
    m2 = maps[1].file_object.map_data()
    if ( (m1.accessor().all() != m2.accessor().all()) or
         (m1.accessor().focus() != m2.accessor().focus()) or
         (m1.accessor().origin() != m2.accessor().origin()) ):
      raise Sorry('The gridding of the two maps is not compatible.')
  else:
  # check if MTZ files have complex arrays and labels
    i_map = 1
    for i in xrange(len(maps)):
      if (maps[i].file_type == 'hkl'):
        labels = get_mtz_labels(maps[i])
        if (len(labels) == 0):
          raise Sorry('%s does not have complex map coefficients' %
                      maps[i].file_name)
        label_phil = getattr(params.input, 'mtz_label_' + str(i_map))
        i_map += 1
        if (label_phil is None):
          raise Sorry('No labels were specified for %s.' % maps[i].file_name)
        elif (label_phil not in labels):
          raise Sorry('%s does not exist in %s' %
                      (label_phil, maps[i].file_name))

  return True

# =============================================================================
# GUI-specific class for running command
from libtbx import runtime_utils
class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    # os.mkdir(self.output_dir)
    # os.chdir(self.output_dir)
    result = run(args=self.args, validated=True)
    return result

# =============================================================================
if (__name__ == "__main__"):
  run(sys.argv[1:])
