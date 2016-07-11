from __future__ import division
import os
from iotbx.phil import parse
from libtbx.utils import Sorry

master_phil_str = """
dispatcher = cctbx.xfel.xtc_process
  .type = str
  .help = Which program to run. cxi.xtc_process is for module only based processing, \
          such as mod_hitfind. cctbx.xfel.xtc_process uses the DIALS back end.
dry_run = False
  .type = bool
  .help = If True, the program will create the trial directory but not submit the job, \
          and will show the command that would have been executed.
experiment = ""
  .type = str
  .help = Experiment name, eg cxid9114
experiment_tag = ""
  .type = str
  .help = User defined tag to describe the set of trials being performed. All database tables will \
          be pre-pended with this string
db {
  host = psdb-user.slac.stanford.edu
    .type = str
    .help = Host name for mysql databse server
  name = ""
    .type = str
    .help = Database name
  user = ""
    .type=str
    .help = Database user name
  password = ""
    .type = str
    .help = Database password. Will be cached as plain text!
}
output_folder = ""
  .type = path
  .help = Processing results will go in this folder
web {
  user = ""
    .type = str
    .help = Username for LCLS run database web service
  password = ""
    .type = str
    .help = Web password. Will be cached in plain text!
  enforce80 = False
    .type = bool
    .help = report only on stream 81, FEE spectrometer
  enforce81 = False
    .type = bool
    .help = report only on stream 81, FEE spectrometer
}
average_raw_data = False
  .type = bool
  .help = If True, don't use any psana corrections (dark, common mode, etc.)

include scope xfel.command_line.cxi_mpi_submit.mp_phil_scope
"""
master_phil_scope = parse(master_phil_str, process_includes=True)

settings_dir = os.path.join(os.path.expanduser('~'), '.cctbx.xfel')
settings_file = os.path.join(settings_dir, 'settings.phil')

def load_cached_settings():
  if os.path.exists(settings_file):
    user_phil = parse(file_name = settings_file)
    return master_phil_scope.fetch(source = user_phil).extract()
  else:
    return master_phil_scope.extract()

def save_cached_settings(params):
  if not os.path.exists(settings_dir):
    os.makedirs(settings_dir)

  working_phil = master_phil_scope.format(python_object = params)
  diff_phil = master_phil_scope.fetch_diff(source = working_phil)

  try:
    f = open(settings_file.encode('utf8'), 'wb')
    f.write(diff_phil.as_str())
    f.close()
  except IOError:
    raise Sorry('Unable to write %s.' % settings_file)
