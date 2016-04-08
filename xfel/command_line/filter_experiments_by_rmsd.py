#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# filter_experiments_by_rmsd.py
#
#  Copyright (C) 2016 Lawrence Berkeley National Laboratory (LBNL)
#
#  Author: Aaron Brewster and David Waterman
#
#  This code is distributed under the X license, a copy of which is
#  included in the root directory of this package.
#
# LIBTBX_SET_DISPATCHER_NAME dev.xfel.filter_experiments_by_rmsd
#
from __future__ import division
from dials.array_family import flex
from scitbx.matrix import col
from matplotlib import pyplot as plt
from libtbx.phil import parse
import libtbx.load_env
import math

help_message = '''
Filter a set of experiments and reflections from a multi-experiment job by overall RMSD
using Tukey's rule of thumb.  I.E, for each experiment, determine the RMSD of the
differences between preditions - observations. Then compute the five number summary of
this set of per-image RMSDs. Then, filter outliers more than iqr_multiplier times the
interquartile range from the third quartile. When x=1.5, this is Tukey's rule.

Example:

  %s combined_experiments.json combined_reflections.pickle
''' % libtbx.env.dispatcher_name

# Create the phil parameters
phil_scope = parse('''
iqr_multiplier = 1.5
  .type = float
  .help = Interquartile multiplier
show_plots = False
  .type = bool
  .help = Show some plots
output {
  filtered_experiments = filtered_experiments.json
    .type = str
    .help = Name of output filtered experiments file
  filtered_reflections = filtered_reflections.pickle
    .type = str
    .help = Name of output filtered reflections file
}
''')

class Script(object):
  ''' Class to parse the command line options. '''

  def __init__(self):
    ''' Set the expected options. '''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # Create the option parser
    usage = "usage: %s combined_experiments.json combined_reflections.pickle" % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      sort_options=True,
      phil=phil_scope,
      read_experiments=True,
      read_reflections=True,
      epilog=help_message)

  def run(self):
    ''' Parse the options. '''
    from dials.util.options import flatten_experiments, flatten_reflections
    from dxtbx.model.experiment.experiment_list import ExperimentList
    from scitbx.math import five_number_summary
    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    self.params = params
    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    assert len(reflections) == 1
    reflections = reflections[0]
    print "Found", len(reflections), "reflections", "and", len(experiments), "experiments"

    difference_vector_norms = (reflections['xyzcal.mm']-reflections['xyzobs.mm.value']).norms()

    data = flex.double()
    counts = flex.double()
    for i in xrange(len(experiments)):
      dvns = difference_vector_norms.select(reflections['id']==i)
      counts.append(len(dvns))
      if len(dvns) == 0:
        data.append(0)
        continue
      rmsd = math.sqrt(flex.sum_sq(dvns)/len(dvns))
      data.append(rmsd)
    data *= 1000
    subset = data.select(counts > 0)
    print len(subset), "experiments with > 0 reflections"

    if params.show_plots:
      h = flex.histogram(subset, n_slots=40)
      fig = plt.figure()
      ax = fig.add_subplot('111')
      ax.plot(h.slot_centers().as_numpy_array(), h.slots().as_numpy_array(), '-')
      plt.title("Histogram of %d image RMSDs"%len(subset))

      fig = plt.figure()
      plt.boxplot(subset, vert=False)
      plt.title("Boxplot of %d image RMSDs"%len(subset))
      plt.show()

    outliers = counts == 0
    min_x, q1_x, med_x, q3_x, max_x = five_number_summary(subset)
    print "Five number summary of RMSDs (microns): min %.1f, q1 %.1f, med %.1f, q3 %.1f, max %.1f"%(min_x, q1_x, med_x, q3_x, max_x)
    iqr_x = q3_x - q1_x
    cut_x = params.iqr_multiplier * iqr_x
    outliers.set_selected(data > q3_x + cut_x, True)
    #outliers.set_selected(col < q1_x - cut_x, True) # Don't throw away the images that are outliers in the 'good' direction!

    filtered_reflections = flex.reflection_table()
    filtered_experiments = ExperimentList()
    for i in xrange(len(experiments)):
      if outliers[i]:
        continue
      refls = reflections.select(reflections['id']==i)
      refls['id'] = flex.int(len(refls), len(filtered_experiments))
      filtered_reflections.extend(refls)
      filtered_experiments.append(experiments[i])

    zeroes = counts == 0
    n_zero = len(counts.select(zeroes))
    print "Removed %d bad experiments and %d experiments with zero reflections, out of %d (%%%.1f)"%(
      len(experiments)-len(filtered_experiments)-n_zero,
      n_zero,
      len(experiments),
      100*((len(experiments)-len(filtered_experiments))/len(experiments)))
    from dxtbx.model.experiment.experiment_list import ExperimentListDumper
    dump = ExperimentListDumper(filtered_experiments)
    dump.as_json(params.output.filtered_experiments)

    filtered_reflections.as_pickle(params.output.filtered_reflections)

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
