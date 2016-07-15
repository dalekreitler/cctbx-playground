from __future__ import division
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.xtc_process
#
try:
  import psana
except ImportError:
  pass # for running at home without psdm build
from xfel.cftbx.detector import cspad_cbf_tbx
from xfel.cxi.cspad_ana import cspad_tbx
import pycbf, os, sys, copy, socket
import libtbx.load_env
from libtbx.utils import Sorry, Usage
from dials.util.options import OptionParser
from libtbx.phil import parse
from dxtbx.imageset import MemImageSet
from dxtbx.datablock import DataBlockFactory
from scitbx.array_family import flex
import numpy as np

xtc_phil_str = '''
  dispatch {
    max_events = None
      .type = int
      .help = If not specified, process all events. Otherwise, only process this many
    process_percent = None
      .type = int(value_min=1, value_max=100)
      .help = Percent of events to process
    estimate_gain_only = False
      .type = bool
      .help = Use to print estimated gain parameters for each event, then exit without attempting \
              further processing.
    hit_finder{
      enable = True
        .type = bool
        .help = Whether to do hitfinding. hit_finder=False: process all images
      minimum_number_of_reflections = 16
        .type = int
        .help = If the number of strong reflections on an image is less than this, and \
                 the hitfinder is enabled, discard this image.
    }
    index = True
      .type = bool
      .help = Attempt to index images
    integrate = True
      .type = bool
      .help = Integrated indexed images. Ignored if index=False
    dump_strong = False
      .type = bool
      .help = Save strongly diffracting images to cbf format
    dump_indexed = True
      .type = bool
      .help = Save indexed images to cbf format
    dump_all = False
      .type = bool
      .help = All frames will be saved to cbf format if set to True
    reindex_strong = False
      .type = bool
      .help = If true, after indexing and refinement, re-index the strong reflections with \
              no outlier rejection
  }
  debug
    .help = Use these flags to track down problematic events that cause unhandled exceptions. \
            Here, a bad event means it caused an unhandled exception, not that the image \
            failed to index. \
            Examples: \
            Process only unprocessed events (excluding bad events): \
              skip_processed_events=True, skip_unprocessed_events=False skip_bad_events=True \
            Process only bad events (for debugging): \
              skip_processed_events=True, skip_unprocessed_events=True skip_bad_events=False \
            Note, due to how MPI works, if an unhandled exception occurrs, some bad events \
            will be marked as bad that were simply in process when the program terminated \
            due to a bad event. Try processing only bad events single process to find the \
            culprit and alert the program authors.
  {
    skip_processed_events = False
      .type = bool
      .help = If True, will look for diagnostic files in the output directory and use \
              them to skip events that had already been processed (succesfully or not)
    skip_unprocessed_events = False
      .type = bool
      .help = If True, will look for diagnostic files in the output directory and use \
              them to skip events that had haven't been processed
    skip_bad_events = False
      .type = bool
      .help = If True, will look for diagnostic files in the output directory and use \
              them to skip events that had caused unhandled exceptions previously
    event_timestamp = None
      .type = str
      .multiple = True
      .help = List of timestamps. If set, will only process the events that match them
  }
  input {
    cfg = None
      .type = str
      .help = Path to psana config file. Genearlly not needed for CBFs. For image pickles, \
              the psana config file should have a mod_image_dict module.
    experiment = None
      .type = str
      .help = Experiment identifier, e.g. cxi84914
    run_num = None
      .type = int
      .help = Run number or run range to process
    address = None
      .type = str
      .help = Detector address, e.g. CxiDs2.0:Cspad.0, or detector alias, e.g. Ds1CsPad
    override_trusted_max = None
      .type = int
      .help = During spot finding, override the saturation value for this data. \
              Overloads will not be integrated, but they can assist with indexing.
    override_trusted_min = None
      .type = int
      .help = During spot finding and indexing, override the minimum pixel value \
              for this data. This does not affect integration.
    use_ffb = False
      .type = bool
      .help = Run on the ffb if possible. Only for active users!
    xtc_dir = None
      .type = str
      .help = Optional path to data directory if it's non-standard. Only needed if xtc \
              streams are not in the standard location for your PSDM installation.
    trial = None
      .type = int
      .help = Optional. Trial number for this run.
    rungroup = None
      .type = int
      .help = Optional. Useful for organizing runs with similar parameters into logical \
              groupings.
 }
  format {
    file_format = *cbf pickle
      .type = choice
      .help = Output file format, 64 tile segmented CBF or image pickle
    pickle {
      out_key = cctbx.xfel.image_dict
        .type = str
        .help = Key name that mod_image_dict uses to put image data in each psana event
    }
    cbf {
      detz_offset = None
        .type = float
        .help = Distance from back of detector rail to sample interaction region (CXI) \
                or actual detector distance (XPP)
      override_energy = None
        .type = float
        .help = If not None, use the input energy for every event instead of the energy \
                from the XTC stream
      override_distance = None
        .type = float
        .help = If not None, use the input distance for every event instead of the distance \
                from the XTC stream
      invalid_pixel_mask = None
        .type = str
        .help = Path to invalid pixel mask, in the dials.generate_mask format. If not set, use the \
                psana computed invalid pixel mask. Regardless, pixels outside of the trusted range \
                for each image will also be masked out. See cxi.make_dials_mask.
      gain_mask_value = None
        .type = float
        .help = If not None, use the gain mask for the run to multiply the low-gain pixels by this number
      common_mode {
        algorithm = default custom
          .type = choice
          .help = Choice of SLAC's common mode correction algorithms. If not specified, use no common \
                  mode correction, only dark pedestal subtraction. Default: use the default common_mode \
                  correction. Custom, see \
                  https://confluence.slac.stanford.edu/display/PSDM/Common+mode+correction+algorithms
        custom_parameterization = None
          .type = ints
          .help = Parameters to control SLAC's common mode correction algorithms. Should be None if \
                  common_mode.algorithm is default or None.  See \
                  https://confluence.slac.stanford.edu/display/PSDM/Common+mode+correction+algorithms
      }
    }
  }
  output {
    output_dir = .
      .type = str
      .help = Directory output files will be placed
    logging_dir = None
      .type = str
      .help = Directory output log files will be placed
    datablock_filename = %s_datablock.json
      .type = str
      .help = The filename for output datablock
    strong_filename = %s_strong.pickle
      .type = str
      .help = The filename for strong reflections from spot finder output.
    indexed_filename = %s_indexed.pickle
      .type = str
      .help = The filename for indexed reflections.
    refined_experiments_filename = %s_refined_experiments.json
      .type = str
      .help = The filename for saving refined experimental models
    integrated_filename = %s_integrated.pickle
      .type = str
      .help = The filename for final integrated reflections.
    profile_filename = None
      .type = str
      .help = The filename for output reflection profile parameters
    integration_pickle = int-%d-%s.pickle
      .type = str
      .help = Filename for cctbx.xfel-style integration pickle files
    reindexedstrong_filename = %s_reindexedstrong.pickle
      .type = str
      .help = The file name for re-indexed strong reflections
  }
  mp {
    method = *mpi sge
      .type = choice
      .help = Muliprocessing method
    mpi {
      method = client_server *striping
        .type = choice
        .help = Method of serving data to child processes in MPI. client_server:    \
                use one process as a server that sends timestamps to each process.  \
                All processes will stay busy at all times at the cost of MPI send/  \
                recieve overhead. striping: each process uses its rank to determine \
                which events to process. Some processes will finish early and go    \
                idle, but no MPI overhead is incurred.
    }
  }
'''

from dials.command_line.stills_process import dials_phil_str

extra_dials_phil_str = '''
  verbosity = 1
   .type = int(value_min=0)
   .help = The verbosity level
  border_mask {
    include scope dials.util.masking.phil_scope
  }
'''

db_logging_phil_str = '''
  experiment_tag = None
    .type = str
    .help = Used if using DB logging. This tag is prepended to each db table name.
  db {
    host = psdb-user.slac.stanford.edu
      .type=str
    name = ""
      .type = str
    user = ""
      .type=str
    password = ""
      .type = str
  }
'''
phil_scope = parse(xtc_phil_str + dials_phil_str + extra_dials_phil_str + db_logging_phil_str, process_includes=True)

# work around for copying dxtbx FormatCBFCspad objects
from xfel.cftbx.detector.cspad_cbf_tbx import cbf_wrapper
def __stupid_but_swig_safe__deepcopy__(self, memo):
  pass
cbf_wrapper.__deepcopy__ = __stupid_but_swig_safe__deepcopy__

from xfel.command_line.xfel_process import Script as DialsProcessScript
class InMemScript(DialsProcessScript):
  """ Script to process XFEL data at LCLS """
  def __init__(self):
    """ Set up the option parser. Arguments come from the command line or a phil file """
    self.usage = """
%s input.experiment=experimentname input.run_num=N input.address=address
 format.file_format=cbf format.cbf.detz_offset=N
%s input.experiment=experimentname input.run_num=N input.address=address
 format.file_format=pickle input.cfg=filename
    """%(libtbx.env.dispatcher_name, libtbx.env.dispatcher_name)
    self.parser = OptionParser(
      usage = self.usage,
      phil = phil_scope)

    self.debug_file_path = None
    self.mpi_log_file_path = None

  def debug_write(self, string):
    debug_file_handle = open(self.debug_file_path, 'a')#, 0)  # 0 for unbuffered
    debug_file_handle.write(string)
    debug_file_handle.close()

  def mpi_log_write(self, string):
    mpi_log_file_handle = open(self.mpi_log_file_path, 'a')#, 0)  # 0 for unbuffered
    mpi_log_file_handle.write(string)
    mpi_log_file_handle.close()

  def run(self):
    """ Process all images assigned to this thread """

    params, options = self.parser.parse_args(
      show_diff_phil=True)

    # Configure the logging
    from dials.util import log
    log.config(params.verbosity)

    # Check inputs
    if params.input.experiment is None or \
       params.input.run_num is None or \
       params.input.address is None:
      raise Usage(self.usage)

    if params.format.file_format == "cbf":
      if params.format.cbf.detz_offset is None:
        raise Usage(self.usage)
    elif params.format.file_format == "pickle":
      if params.input.cfg is None:
        raise Usage(self.usage)
    else:
      raise Usage(self.usage)

    if not os.path.exists(params.output.output_dir):
      raise Sorry("Output path not found:" + params.output.output_dir)

    # The convention is to put %s in the phil parameter to add a time stamp to
    # each output datafile. Save the initial templates here.
    self.strong_filename_template              = params.output.strong_filename
    self.indexed_filename_template             = params.output.indexed_filename
    self.refined_experiments_filename_template = params.output.refined_experiments_filename
    self.integrated_filename_template          = params.output.integrated_filename
    self.reindexedstrong_filename_template     = params.output.reindexedstrong_filename

    # Don't allow the strong reflections to be written unless there are enough to
    # process
    params.output.strong_filename = None

    # Save the paramters
    self.params_cache = copy.deepcopy(params)
    self.options = options

    if params.mp.method == "mpi":
      from mpi4py import MPI
      comm = MPI.COMM_WORLD
      rank = comm.Get_rank() # each process in MPI has a unique id, 0-indexed
      size = comm.Get_size() # size: number of processes running in this job
    elif params.mp.method == "sge" and \
        'SGE_TASK_ID'    in os.environ and \
        'SGE_TASK_FIRST' in os.environ and \
        'SGE_TASK_LAST'  in os.environ:
      if 'SGE_STEP_SIZE' in os.environ:
        assert int(os.environ['SGE_STEP_SIZE']) == 1
      if os.environ['SGE_TASK_ID'] == 'undefined' or os.environ['SGE_TASK_ID'] == 'undefined' or os.environ['SGE_TASK_ID'] == 'undefined':
        rank = 0
        size = 1
      else:
        rank = int(os.environ['SGE_TASK_ID']) - int(os.environ['SGE_TASK_FIRST'])
        size = int(os.environ['SGE_TASK_LAST']) - int(os.environ['SGE_TASK_FIRST']) + 1
    else:
      rank = 0
      size = 1

    if params.output.logging_dir is not None:
      log_path = os.path.join(params.output.logging_dir, "log_rank%04d.out"%rank)
      error_path = os.path.join(params.output.logging_dir, "error_rank%04d.out"%rank)
      print "Redirecting stdout to %s"%log_path
      print "Redirecting stderr to %s"%error_path
      sys.stdout = open(log_path,'a', buffering=0)
      sys.stderr = open(error_path,'a',buffering=0)
      print "Should be redirected now"

    debug_dir = os.path.join(params.output.output_dir, "debug")
    if not os.path.exists(debug_dir):
      try:
        os.makedirs(debug_dir)
      except OSError, e:
        pass # due to multiprocessing, makedirs can sometimes fail
    assert os.path.exists(debug_dir)

    if params.debug.skip_processed_events or params.debug.skip_processed_events or params.debug.skip_bad_events:
      print "Reading debug files..."
      self.known_events = {}
      for filename in os.listdir(debug_dir):
        # format: hostname,timestamp,status
        for line in open(os.path.join(debug_dir, filename)):
          vals = line.strip().split(',')
          if len(vals) == 2:
            self.known_events[vals[1]] = "unknown"
          elif len(vals) == 3:
            self.known_events[vals[1]] = vals[2]

    self.debug_file_path = os.path.join(debug_dir, "debug_%d.txt"%rank)
    write_newline = os.path.exists(self.debug_file_path)
    if write_newline: # needed if the there was a crash
      self.debug_write("\n")

    if params.mp.method != 'mpi' or params.mp.mpi.method == 'client_server':
      if rank == 0:
        self.mpi_log_file_path = os.path.join(debug_dir, "mpilog.out")
        write_newline = os.path.exists(self.mpi_log_file_path)
        if write_newline: # needed if the there was a crash
          self.mpi_log_write("\n")

    # set up psana
    if params.input.cfg is not None:
      psana.setConfigFile(params.input.cfg)
    dataset_name = "exp=%s:run=%s:idx"%(params.input.experiment,params.input.run_num)
    if params.input.xtc_dir is not None:
      if params.input.use_ffb:
        raise Sorry("Cannot specify the xtc_dir and use SLAC's ffb system")
      dataset_name += ":dir=%s"%params.input.xtc_dir
    elif params.input.use_ffb:
      # as ffb is only at SLAC, ok to hardcode /reg/d here
      dataset_name += ":dir=/reg/d/ffb/%s/%s/xtc"%(params.input.experiment[0:3],params.input.experiment)
    ds = psana.DataSource(dataset_name)

    if params.format.file_format == "cbf":
      self.psana_det = psana.Detector(params.input.address, ds.env())

    # set this to sys.maxint to analyze all events
    if params.dispatch.max_events is None:
      max_events = sys.maxint
    else:
      max_events = params.dispatch.max_events

    for run in ds.runs():
      if params.format.file_format == "cbf":
        # load a header only cspad cbf from the slac metrology
        self.base_dxtbx = cspad_cbf_tbx.env_dxtbx_from_slac_metrology(run, params.input.address)
        if self.base_dxtbx is None:
          raise Sorry("Couldn't load calibration file for run %d"%run.run())

        if params.format.file_format == 'cbf':
          if params.format.cbf.common_mode.algorithm == "custom":
            self.common_mode = params.format.cbf.common_mode.custom_parameterization
            assert self.common_mode is not None
          else:
            self.common_mode = params.format.cbf.common_mode.algorithm # could be None or default

        if params.format.cbf.invalid_pixel_mask is not None:
          from libtbx import easy_pickle
          self.psana_mask = easy_pickle.load(params.format.cbf.invalid_pixel_mask)
          assert len(self.psana_mask) == 64
        else:
          psana_mask = self.psana_det.mask(run,calib=True,status=True,edges=True,central=True,unbond=True,unbondnbrs=True)
          if psana_mask.dtype == np.bool:
            psana_mask = flex.bool(psana_mask)
          else:
            psana_mask = flex.bool(psana_mask == 1)
          assert psana_mask.focus() == (32, 185, 388)
          self.psana_mask = []
          for i in xrange(32):
            self.psana_mask.append(psana_mask[i:i+1,:,:194])
            self.psana_mask[-1].reshape(flex.grid(185,194))
            self.psana_mask.append(psana_mask[i:i+1,:,194:])
            self.psana_mask[-1].reshape(flex.grid(185,194))

      # list of all events
      times = run.times()
      nevents = min(len(times),max_events)
      times = times[:nevents]
      if params.dispatch.process_percent is not None:
        import fractions
        percent = params.dispatch.process_percent / 100
        f = fractions.Fraction(percent).limit_denominator(100)
        times = [times[i] for i in xrange(len(times)) if i % f.denominator < f.numerator]
        print "Dividing %d of %d events (%4.1f%%) between all processes"%(len(times), nevents, 100*len(times)/nevents)
        nevents = len(times)
      else:
        print "Dividing %d events between all processes" % nevents
      if params.mp.method == "mpi" and params.mp.mpi.method == 'client_server' and size > 2:
        print "Using MPI client server"
        # use a client/server approach to be sure every process is busy as much as possible
        # only do this if there are more than 2 processes, as one process will be a server
        try:
          if rank == 0:
            # server process
            self.mpi_log_write("MPI START\n")
            for t in times:
              # a client process will indicate it's ready by sending its rank
              self.mpi_log_write("Getting next available process\n")
              rankreq = comm.recv(source=MPI.ANY_SOURCE)
              ts = cspad_tbx.evt_timestamp((t.seconds(),t.nanoseconds()/1e6))
              self.mpi_log_write("Process %s is ready, sending ts %s\n"%(rankreq, ts))
              comm.send(t,dest=rankreq)
            # send a stop command to each process
            self.mpi_log_write("MPI DONE, sending stops\n")
            for rankreq in range(size-1):
              self.mpi_log_write("Getting next available process\n")
              rankreq = comm.recv(source=MPI.ANY_SOURCE)
              self.mpi_log_write("Sending stop to %d\n"%rankreq)
              comm.send('endrun',dest=rankreq)
          else:
            # client process
            while True:
              # inform the server this process is ready for an event
              comm.send(rank,dest=0)
              evttime = comm.recv(source=0)
              if evttime == 'endrun': break
              self.process_event(run, evttime)
        except Exception, e:
          print "Error caught in main loop"
          print str(e)
      else:
        # chop the list into pieces, depending on rank.  This assigns each process
        # events such that the get every Nth event where N is the number of processes
        print "Striping events"
        mytimes = [times[i] for i in xrange(len(times)) if (i+rank)%size == 0]

        for i in xrange(len(mytimes)):
          self.process_event(run, mytimes[i])

      run.end()
    ds.end()

  def process_event(self, run, timestamp):
    """
    Process a single event from a run
    @param run psana run object
    @param timestamp psana timestamp object
    """
    ts = cspad_tbx.evt_timestamp((timestamp.seconds(),timestamp.nanoseconds()/1e6))
    if ts is None:
      print "No timestamp, skipping shot"
      return

    if len(self.params_cache.debug.event_timestamp) > 0 and ts not in self.params_cache.debug.event_timestamp:
      return

    if self.params_cache.debug.skip_processed_events or self.params_cache.debug.skip_unprocessed_events or self.params_cache.debug.skip_bad_events:
      if ts in self.known_events:
        if self.known_events[ts] == "unknown":
          if self.params_cache.debug.skip_bad_events and self.known_events[ts] == "unknown":
            print "Skipping event %s: possibly caused an unknown exception previously"%ts
            return
        elif self.params_cache.debug.skip_processed_events:
          print "Skipping event %s: processed successfully previously"%ts
          return
      else:
        if self.params_cache.debug.skip_unprocessed_events:
          print "Skipping event %s: not processed previously"%ts
          return

    self.debug_write("%s,%s"%(socket.gethostname(), ts))

    evt = run.event(timestamp)
    if evt.get("skip_event") or "skip_event" in [key.key() for key in evt.keys()]:
      print "Skipping event",ts
      self.debug_write(",psana_skip\n")
      return

    print "Accepted", ts
    self.params = copy.deepcopy(self.params_cache)

    # the data needs to have already been processed and put into the event by psana
    if self.params.format.file_format == 'cbf':
      # get numpy array, 32x185x388
      data = cspad_cbf_tbx.get_psana_corrected_data(self.psana_det, evt, use_default=False, dark=True,
                                                    common_mode=self.common_mode,
                                                    apply_gain_mask=self.params.format.cbf.gain_mask_value is not None,
                                                    gain_mask_value=self.params.format.cbf.gain_mask_value,
                                                    per_pixel_gain=False)
      if data is None:
        print "No data"
        self.debug_write(",no_data\n")
        return

      if self.params.format.cbf.override_distance is None:
        distance = cspad_tbx.env_distance(self.params.input.address, run.env(), self.params.format.cbf.detz_offset)
        if distance is None:
          print "No distance, skipping shot"
          self.debug_write(",no_distance\n")
          return
      else:
        distance = self.params.format.cbf.override_distance

      if self.params.format.cbf.override_energy is None:
        wavelength = cspad_tbx.evt_wavelength(evt)
        if wavelength is None:
          print "No wavelength, skipping shot"
          self.debug_write(",no_wavelength\n")
          return
      else:
        wavelength = 12398.4187/self.params.format.cbf.override_energy

    if self.params.format.file_format == 'pickle':
      image_dict = evt.get(self.params.format.pickle.out_key)
      data = image_dict['DATA']

    timestamp = t = ts
    s = t[0:4] + t[5:7] + t[8:10] + t[11:13] + t[14:16] + t[17:19] + t[20:23]
    print "Processing shot", s

    if self.params.format.file_format == 'cbf':
      # stitch together the header, data and metadata into the final dxtbx format object
      cspad_img = cspad_cbf_tbx.format_object_from_data(self.base_dxtbx, data, distance, wavelength, timestamp, self.params.input.address)
    elif self.params.format.file_format == 'pickle':
      from dxtbx.format.FormatPYunspecifiedStill import FormatPYunspecifiedStillInMemory
      cspad_img = FormatPYunspecifiedStillInMemory(image_dict)

    cspad_img.timestamp = s

    if self.params.dispatch.dump_all:
      self.save_image(cspad_img, self.params, os.path.join(self.params.output.output_dir, "shot-" + s))

    self.cache_ranges(cspad_img, self.params)

    imgset = MemImageSet([cspad_img])
    if self.params.dispatch.estimate_gain_only:
      from dials.command_line.estimate_gain import estimate_gain
      estimate_gain(imgset)
      return

    datablock = DataBlockFactory.from_imageset(imgset)[0]

    # before calling DIALS for processing, set output paths according to the templates
    if self.indexed_filename_template is not None and "%s" in self.indexed_filename_template:
      self.params.output.indexed_filename = os.path.join(self.params.output.output_dir, self.indexed_filename_template%("idx-" + s))
    if "%s" in self.refined_experiments_filename_template:
      self.params.output.refined_experiments_filename = os.path.join(self.params.output.output_dir, self.refined_experiments_filename_template%("idx-" + s))
    if "%s" in self.integrated_filename_template:
      self.params.output.integrated_filename = os.path.join(self.params.output.output_dir, self.integrated_filename_template%("idx-" + s))
    if "%s" in self.reindexedstrong_filename_template:
      self.params.output.reindexedstrong_filename = os.path.join(self.params.output.output_dir, self.reindexedstrong_filename_template%("idx-" + s))

    # Load a dials mask from the trusted range and psana mask
    from dials.util.masking import MaskGenerator
    generator = MaskGenerator(self.params.border_mask)
    mask = generator.generate(imgset)
    if self.params.format.file_format == "cbf":
      mask = tuple([a&b for a, b in zip(mask,self.psana_mask)])

    self.params.spotfinder.lookup.mask = mask
    self.params.integration.lookup.mask = mask

    try:
      observed = self.find_spots(datablock)
    except Exception, e:
      import traceback; traceback.print_exc()
      print str(e), "event", timestamp
      self.debug_write(",spotfinding_exception\n")
      return

    print "Found %d bright spots"%len(observed)

    if self.params.dispatch.hit_finder.enable and len(observed) < self.params.dispatch.hit_finder.minimum_number_of_reflections:
      print "Not enough spots to index"
      self.debug_write(",not_enough_spots_%d\n"%len(observed))
      self.log_frame(None, None, run.run(), len(observed), timestamp)
      return

    self.restore_ranges(cspad_img, self.params)

    # save cbf file
    if self.params.dispatch.dump_strong:
      self.save_image(cspad_img, self.params, os.path.join(self.params.output.output_dir, "hit-" + s))

      # save strong reflections.  self.find_spots() would have done this, but we only
      # want to save data if it is enough to try and index it
      if self.strong_filename_template:
        if "%s" in self.strong_filename_template:
          strong_filename = self.strong_filename_template%("hit-" + s)
        else:
          strong_filename = self.strong_filename_template
        strong_filename = os.path.join(self.params.output.output_dir, strong_filename)

        from dials.util.command_line import Command
        Command.start('Saving {0} reflections to {1}'.format(
            len(observed), os.path.basename(strong_filename)))
        observed.as_pickle(strong_filename)
        Command.end('Saved {0} observed to {1}'.format(
            len(observed), os.path.basename(strong_filename)))

    if not self.params.dispatch.index:
      self.debug_write(",strong_shot_%d\n"%len(observed))
      self.log_frame(None, None, run.run(), len(observed), timestamp)
      return

    # index and refine
    try:
      experiments, indexed = self.index(datablock, observed)
    except Exception, e:
      import traceback; traceback.print_exc()
      print str(e), "event", timestamp
      self.debug_write(",indexing_failed_%d\n"%len(observed))
      self.log_frame(None, None, run.run(), len(observed), timestamp)
      return

    if self.params.dispatch.dump_indexed:
      self.save_image(cspad_img, self.params, os.path.join(self.params.output.output_dir, "idx-" + s))

    try:
      experiments = self.refine(experiments, indexed)
    except Exception, e:
      import traceback; traceback.print_exc()
      print str(e), "event", timestamp
      self.debug_write(",refine_failed_%d\n"%len(indexed))
      self.log_frame(None, None, run.run(), len(observed), timestamp)
      return

    if self.params.dispatch.reindex_strong:
      try:
        self.reindex_strong(experiments, observed)
      except Exception, e:
        import traceback; traceback.print_exc()
        print str(e), "event", timestamp
        self.debug_write(",reindexstrong_failed_%d\n"%len(indexed))
        self.log_frame(None, None, run.run(), len(observed), timestamp)
        return

    if not self.params.dispatch.integrate:
      self.debug_write(",index_ok_%d\n"%len(indexed))
      self.log_frame(None, None, run.run(), len(observed), timestamp)
      return

    # integrate
    try:
      integrated = self.integrate(experiments, indexed)
    except Exception, e:
      import traceback; traceback.print_exc()
      print str(e), "event", timestamp
      self.debug_write(",integrate_failed_%d\n"%len(indexed))
      self.log_frame(None, None, run.run(), len(observed), timestamp)
      return

    self.log_frame(experiments, integrated, run.run(), len(observed), timestamp)
    self.debug_write(",integrate_ok_%d\n"%len(integrated))

  def log_frame(self, experiments, reflections, run, n_strong, timestamp = None):
    if self.params.experiment_tag is None:
      return
    try:
      from xfel.ui.db.dxtbx_db import log_frame
      log_frame(experiments, reflections, self.params, run, n_strong, timestamp)
    except Exception, e:
      import traceback; traceback.print_exc()
      print str(e), "event", timestamp
      self.debug_write(",db_logging_failed_%d\n" % len(integrated))

  def save_image(self, image, params, root_path):
    """ Save an image, in either cbf or pickle format.
    @param image dxtbx format object
    @param params phil scope object
    @param root_path output file path without extension
    """

    if params.format.file_format == 'cbf':
      dest_path = root_path + ".cbf"
    elif params.format.file_format == 'pickle':
      dest_path = root_path + ".pickle"

    try:
      if params.format.file_format == 'cbf':
        image._cbf_handle.write_widefile(dest_path, pycbf.CBF,\
          pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K, 0)
      elif params.format.file_format == 'pickle':
        from libtbx import easy_pickle
        easy_pickle.dump(dest_path, image._image_file)
    except Exception:
      print "Warning, couldn't save image:", dest_path

  def cache_ranges(self, cspad_img, params):
    """ Save the current trusted ranges, and replace them with the given overrides, if present.
    @param cspad_image dxtbx format object
    @param params phil scope
    """
    if params.input.override_trusted_max is None and params.input.override_trusted_min is None:
      return

    detector = cspad_img.get_detector()
    self.cached_ranges = []
    for panel in detector:
      new_range = cached_range = panel.get_trusted_range()
      self.cached_ranges.append(cached_range)
      if params.input.override_trusted_max is not None:
        new_range = new_range[0], params.input.override_trusted_max
      if params.input.override_trusted_min is not None:
        new_range = params.input.override_trusted_min, new_range[1]

      panel.set_trusted_range(new_range)

  def restore_ranges(self, cspad_img, params):
    """ Restore the previously cached trusted ranges, if present.
    @param cspad_image dxtbx format object
    @param params phil scope
    """
    if params.input.override_trusted_max is None and params.input.override_trusted_min is None:
      return

    detector = cspad_img.get_detector()
    for cached_range, panel in zip(self.cached_ranges, detector):
      panel.set_trusted_range(cached_range)

  def reindex_strong(self, experiments, strong):
    print "Reindexing strong reflections using refined experimental models and no outlier rejection..."
    from dials.algorithms.indexing.stills_indexer import stills_indexer_known_orientation
    indexer = stills_indexer_known_orientation(strong, experiments.imagesets(), self.params, experiments.crystals())
    indexed_reflections = indexer.reflections.select(indexer.indexed_reflections)

    print "Indexed %d strong reflections out of %d"%(len(indexed_reflections), len(strong))
    self.save_reflections(indexed_reflections, self.params.output.reindexedstrong_filename)

if __name__ == "__main__":
  from dials.util import halraiser
  try:
    script = InMemScript()
    script.run()
  except Exception as e:
    halraiser(e)
