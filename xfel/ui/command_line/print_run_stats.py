from __future__ import division

'''
Author      : Young, I.D.
Created     : 07/14/2016
Last Changed: 07/14/2016
Description : XFEL UI plot real-time run stats
'''

from libtbx.phil import parse
from libtbx.utils import Sorry
from xfel.ui.db.xfel_db import xfel_db_application
from xfel.ui.db.stats import HitrateStats
import sys
from xfel.ui.command_line.plot_run_stats import phil_scope

def run(args):
  user_phil = []
  for arg in args:
    try:
      user_phil.append(parse(arg))
    except Exception, e:
      raise Sorry("Unrecognized argument %s"%arg)
  params = phil_scope.fetch(sources=user_phil).extract()
  print "Printing results for trial", params.trial, "using a hit cutoff of", params.hit_cutoff, "reflections"
  print
  print " Run   N Hits   (%) N Indexed   (%) N High qual   (%)  %HQR   N Frames"

  hit_total = 0
  indexed_total = 0
  high_quality_total = 0
  overall_total = 0

  app = xfel_db_application(params)
  for run_no in params.run:
    try:
      timestamps, n_strong, average_i_sigi_low, average_i_sigi_high = HitrateStats(app, run_no, params.trial, params.rungroup, params.d_min)()
    except Exception, e:
      print "Coulnd't get run", run_no
      continue
    n_hit = (n_strong >= params.hit_cutoff).count(True)
    n_indexed = (average_i_sigi_low > 0).count(True)
    n_total = len(timestamps)
    n_high_quality = (average_i_sigi_high > 0).count(True)
    try:
      print "% 4d  % 7d % 5.1f   % 7d % 5.1f     % 7d % 5.1f % 5.1f    % 7d " % (run_no, n_hit, 100*n_hit/n_total, n_indexed, 100*n_indexed/n_total, n_high_quality, 100*n_high_quality/n_total, 100*n_high_quality/n_indexed, n_total)
    except ZeroDivisionError:
      print "% 4d  % 7d % 5.1f   % 7d % 5.1f     % 7d % 5.1f % 5.1f    % 7d " % (run_no, n_hit, 0, n_indexed, 0, n_high_quality, 0, 0, n_total)

    hit_total += n_hit
    indexed_total += n_indexed
    high_quality_total += n_high_quality
    overall_total += n_total

  if len(params.run) > 1:
    print "-" * 80
    try:
      print "Total % 7d % 5.1f   % 7d % 5.1f     % 7d % 5.1f % 5.1f    % 7d " % (hit_total, 100*hit_total/overall_total, indexed_total, 100*indexed_total/overall_total, high_quality_total, 100*high_quality_total/overall_total, 100*high_quality_total/indexed_total, overall_total)
    except ZeroDivisionError:
      if overall_total == 0:
        print "Total % 7d % 5.1f   % 7d % 5.1f     % 7d % 5.1f % 5.1f    % 7d " % (hit_total, 0, indexed_total, 0, 0, 0, 0, overall_total)
      else:
        print "Total % 7d % 5.1f   % 7d % 5.1f     % 7d % 5.1f % 5.1f    % 7d " % (hit_total, 100*hit_total/overall_total, indexed_total, 100*indexed_total/overall_total, high_quality_total, 100*high_quality_total/overall_total, 0, overall_total)

if __name__ == "__main__":
  run(sys.argv[1:])
