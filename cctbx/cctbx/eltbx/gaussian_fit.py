from cctbx.eltbx import xray_scattering
from cctbx.array_family import flex
from scitbx import lbfgs
from scitbx.python_utils import easy_pickle
from scitbx.python_utils.misc import adopt_init_args, user_plus_sys_time
import math
import sys, os

international_tables_sampling_points_stol = flex.double(
  [0.00, 0.01, 0.02, 0.03, 0.04, 0.05,
   0.06, 0.07, 0.08, 0.09, 0.10,
   0.11, 0.12, 0.13, 0.14, 0.15,
   0.16, 0.17, 0.18, 0.19, 0.20,
   0.22, 0.24, 0.25, 0.26, 0.28, 0.30,
   0.32, 0.34, 0.35, 0.36, 0.38, 0.40,
   0.42, 0.44, 0.45, 0.46, 0.48, 0.50,
   0.55, 0.60, 0.65, 0.70, 0.80, 0.90, 1.00,
   1.10, 1.20, 1.30, 1.40, 1.50,
   1.60, 1.70, 1.80, 1.90, 2.00,
   2.50, 3.00, 3.50, 4.00, 5.00, 6.00])

international_tables_sampling_points_d_star_sq = \
  flex.pow2(international_tables_sampling_points_stol * 2)

international_tables_sampled_value_sigmas = flex.double(
  [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005,
   0.0005, 0.0005, 0.0005, 0.0005, 0.0005,
   0.0005, 0.0005, 0.0005, 0.0005, 0.0005,
   0.0005, 0.0005, 0.0005, 0.0005, 0.0005,
   0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005,
   0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005,
   0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005,
   0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005,
   0.0005, 0.0005, 0.0005, 0.0005, 0.0005,
   0.0005, 0.0005, 0.0005, 0.0005, 0.0005,
   0.0008, 0.0008, 0.0008, 0.0008, 0.0008, 0.0008])

class international_tables_sampling_points_and_value_sigmas_up_to:

  def __init__(self, d_min, eps=1.e-6):
    assert d_min > 0
    self.selection = international_tables_sampling_points_d_star_sq \
                     < (1/d_min)**2 + eps
    self.d_star_sq = international_tables_sampling_points_d_star_sq.select(
      self.selection)
    self.sigmas = international_tables_sampled_value_sigmas.select(
      self.selection)

class minimize:

  def __init__(self, diff_gaussian, target_power, d_star_sq,
                     weights=None,
                     b_min=-1,
                     enforce_positive_b=0001,
                     lbfgs_termination_params=None,
                     lbfgs_core_params=lbfgs.core_parameters(m=7)):
    adopt_init_args(self, locals())
    assert target_power in [2,4]
    self.n = diff_gaussian.n_ab() * 2
    if (weights is None):
      self.weights = flex.double(d_star_sq.size(), 1)
    self.x = flex.double(self.n, 0)
    self.first_target_value = None
    self.minimizer = lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      core_params=lbfgs_core_params)
    self.apply_shifts()
    self.compute_target(compute_gradients=00000)
    self.final_target_value = self.f
    self.final_diff_gaussian = self.diff_gaussian_shifted

  def apply_shifts(self):
    if (not self.enforce_positive_b):
      shifts = self.x
    else:
      n_ab = self.diff_gaussian.n_ab()
      b = flex.double(self.diff_gaussian.b())
      assert flex.min(b) >= 0
      sqrt_b = flex.sqrt(b)
      shifts_sqrt_b = self.x[n_ab:]
      sqrt_b_shifted = sqrt_b + shifts_sqrt_b
      b_shifted = flex.pow2(sqrt_b_shifted)
      shifts_b = b_shifted - b
      shifts = self.x[:n_ab]
      shifts.append(shifts_b)
      assert min(self.diff_gaussian.b()) >= 0
    self.diff_gaussian_shifted = self.diff_gaussian.apply_shifts(
      shifts, self.b_min)

  def compute_target(self, compute_gradients):
    dg = self.diff_gaussian_shifted
    tt = dg.target_terms_at_points(self.d_star_sq)
    self.f = flex.pow2(tt)
    if (self.target_power == 4):
      self.f = flex.pow2(self.f)
    self.f = flex.sum(self.weights * self.f)
    if (compute_gradients):
      self.g = dg.sum_of_gradients_at_points(
        self.target_power, self.d_star_sq, self.weights, tt, 00000)
      if (self.enforce_positive_b):
        n_ab = self.diff_gaussian.n_ab()
        b = flex.double(self.diff_gaussian.b())
        assert flex.min(b) >= 0
        sqrt_b = flex.sqrt(b)
        shifts_sqrt_b = self.x[n_ab:]
        d_b_d_shift = (sqrt_b + shifts_sqrt_b) * 2
        g_shifts_sqrt_b = self.g[n_ab:] * d_b_d_shift
        self.g = self.g[:n_ab]
        self.g.append(g_shifts_sqrt_b)
    else:
      self.g = None

  def __call__(self):
    if (self.first_target_value is None):
      assert self.x.all_eq(0)
      self.diff_gaussian_shifted = self.diff_gaussian
    else:
      self.apply_shifts()
    self.compute_target(compute_gradients=0001)
    if (self.first_target_value is None):
      self.first_target_value = self.f
    return self.x, self.f, self.g

class d_interval:

  def __init__(self, d_max, d_min):
    assert 0 < d_min < d_max
    self.d_max = d_max
    self.d_min = d_min
    self.d_star_min = 1./d_max
    self.d_star_max = 1./d_min
    self.d_star_range = self.d_star_max - self.d_star_min

  def d_star(self, f_range):
    return self.d_star_min + f_range * self.d_star_range

  def d_star_sq(self, f_range):
    return self.d_star(f_range)**2

  def stol_sq(self, f_range):
    return self.d_star_sq(f_range) / 4

  def d(self, f_range):
    return 1./self.d_star(f_range)

def make_start_gaussian(reference_gaussian,
                        existing_gaussian,
                        start_interval,
                        i_start_fraction,
                        f_range):
  assert reference_gaussian.n_ab() >= existing_gaussian.n_ab()
  stol_sq = start_interval.stol_sq(f_range)
  f0_reference = reference_gaussian.at_stol_sq(0)
  fs_reference = reference_gaussian.at_stol_sq(stol_sq)
  f0_existing = existing_gaussian.at_stol_sq(0)
  fs_existing = existing_gaussian.at_stol_sq(stol_sq)
  n_terms = existing_gaussian.n_ab() + 1
  if (n_terms == 1):
    a = flex.double([f0_reference])
    b = flex.double()
    fs_part = fs_reference
  else:
    scale_old = (1-(1.+i_start_fraction)/n_terms)
    a = flex.double(existing_gaussian.a()) * scale_old
    a.append(f0_reference - flex.sum(a))
    b = flex.double(existing_gaussian.b())
    fs_part = fs_reference - fs_existing * scale_old
  addl_b = 0
  if (a[-1] != 0):
    r = fs_part / a[-1]
    if (0 < r <= 1):
      addl_b = -math.log(r) / stol_sq
  b.append(addl_b)
  if (addl_b != 0):
    assert abs(a[-1] * math.exp(-b[-1] * stol_sq) - fs_part) < 1.e-6
  result = xray_scattering.difference_gaussian(
    reference_gaussian,
    xray_scattering.gaussian(iter(a), iter(b)))
  if (addl_b != 0):
    assert abs(result.at_stol_sq(0) - f0_reference) < 1.e-4
  if (n_terms == 1):
    assert abs(result.at_stol_sq(stol_sq) - fs_reference) < 1.e-4
  return result

def get_significant_relative_errors(diff_gaussian, d_star_sq, sigmas):
  assert d_star_sq.size() == sigmas.size()
  reference_values = diff_gaussian.reference_gaussian().at_d_star_sq(d_star_sq)
  diffs = diff_gaussian.target_terms_at_points(d_star_sq)
  results = flex.double()
  for i,diff in diffs.items():
    sigma = sigmas[i]
    result = max(0, abs(diff)-sigma)
    if (result > 0):
      assert abs(reference_values[i]) != 0
      result /= abs(reference_values[i])
    results.append(result)
  return results

class find_d_min:

  def __init__(self, diff_gaussian,
                     target_power,
                     max_max_error,
                     start_interval,
                     min_interval):
    interval = [start_interval.d_max, start_interval.d_min]
    good_d_min = None
    while 1:
      self.d_min = (interval[0] + interval[1]) / 2.
      points = international_tables_sampling_points_and_value_sigmas_up_to(
        d_min=self.d_min)
      self.min = minimize(
        diff_gaussian=diff_gaussian,
        target_power=target_power,
        d_star_sq=points.d_star_sq)
      self.max_error = flex.max(get_significant_relative_errors(
        diff_gaussian=self.min.final_diff_gaussian,
        d_star_sq=points.d_star_sq,
        sigmas=points.sigmas))
      if (    self.max_error <= max_max_error
          and min(self.min.final_diff_gaussian.b()) > self.min.b_min):
        if (interval[0] == start_interval.d_min):
          break
        if (interval[0] - interval[1] <= min_interval):
          if (interval[1] - start_interval.d_min <= min_interval):
            interval[0] = start_interval.d_min
            interval[1] = start_interval.d_min
          else:
            break
        else:
          interval[0] = self.d_min
        diff_gaussian = self.min.final_diff_gaussian
        good_d_min = self.d_min
        good_min = self.min
        good_max_error = self.max_error
      else:
        if (interval[0] - interval[1] <= min_interval):
          if (good_d_min is not None):
            self.d_min = good_d_min
            self.min = good_min
            self.max_error = good_max_error
          else:
            self.d_min = None
            self.min = None
            self.max_error = None
          break
        interval[1] = self.d_min

class find_d_min_multi:

  def __init__(self, reference_gaussian,
                     existing_gaussian,
                     target_powers,
                     start_interval,
                     n_trial_points_factor,
                     n_start_fractions,
                     max_max_error,
                     min_interval):
    n_terms = existing_gaussian.n_ab() + 1
    n_trial_points = n_terms * n_trial_points_factor
    best_fit = None
    for i_trial_point in xrange(n_trial_points+1):
      f_range = i_trial_point / float(n_trial_points)
      d_min = start_interval.d(f_range)
      for i_start_fraction in xrange(min(n_start_fractions,n_terms)):
        diff_gaussian = make_start_gaussian(
          reference_gaussian=reference_gaussian,
          existing_gaussian=existing_gaussian,
          start_interval=start_interval,
          i_start_fraction=i_start_fraction,
          f_range=f_range/2)
        for target_power in target_powers:
          fit = find_d_min(
            diff_gaussian=diff_gaussian,
            target_power=target_power,
            max_max_error=max_max_error,
            start_interval=start_interval,
            min_interval=min_interval)
          if (fit.d_min is not None):
            if (best_fit is None
                or best_fit.d_min > fit.d_min
                or best_fit.d_min == fit.d_min
                  and best_fit.max_error > fit.max_error):
              best_fit = fit
    if (best_fit is None):
      self.d_min = None
      self.min = None
      self.max_error = None
    else:
      self.d_min = best_fit.d_min
      self.min = best_fit.min
      self.max_error = best_fit.max_error

def show_fit_summary(source, label, gaussian, d_min, e, e_other=None):
  n_terms = str(gaussian.n_ab())
  if (gaussian.c() != 0): n_terms += "+c"
  n_terms += ","
  print "%24s: %s n_terms=%-4s d_min=%.2f, e=%.4f" % (
    source, label, n_terms, d_min, e),
  if (e_other is not None and e_other > e):
    print "Better",
  print

def show_literature_fits(label, n_terms, reference_gaussian,
                         sampling_points, d_min, e_other=None):
  for lib in [xray_scattering.it1992,
              xray_scattering.two_gaussian_agarwal_isaacs,
              xray_scattering.two_gaussian_agarwal_1978,
              xray_scattering.one_gaussian_agarwal_1978]:
    if (lib == xray_scattering.it1992):
      lib_gaussian = xray_scattering.it1992(label, 1).fetch()
      lib_source = "IT1992"
    elif (lib.table.has_key(label)):
      lib_gaussian = lib.table[label]
      lib_source = lib.source_short
    else:
      lib_gaussian = None
    if (lib_gaussian is not None and lib_gaussian.n_ab() == n_terms):
      diff_gaussian = xray_scattering.difference_gaussian(
        reference_gaussian, lib_gaussian)
      e = flex.max(get_significant_relative_errors(
        diff_gaussian=diff_gaussian,
        d_star_sq=sampling_points.d_star_sq,
        sigmas=sampling_points.sigmas))
      show_fit_summary(lib_source, label, lib_gaussian, d_min, e, e_other)

def write_plot(f, gaussian, d_star_sq):
  for dss in d_star_sq:
    d_star = math.sqrt(dss)
    print >> f, d_star, gaussian.at_d_star_sq(dss)
  print >> f, "&"

def write_plots(plots_dir, label, gaussians, d_star_sq):
  label = label.replace("'", "prime")
  file_name = os.path.join(plots_dir, label+".xy")
  f = open(file_name, "w")
  for gaussian in gaussians:
    write_plot(f, gaussian, d_star_sq)
  f.close()

class fit_parameters:

  def __init__(self, max_n_terms=5,
                     target_powers=[2,4],
                     max_max_error=0.01,
                     start_interval=d_interval(d_max=15, d_min=1/12.),
                     find_d_min_min_interval=0.01,
                     n_trial_points_factor=5,
                     n_start_fractions=3):
    adopt_init_args(self, locals())

def incremental_fits(label, reference_gaussian, params=None, plots_dir=None,
                     verbose=0):
  if (params is None): params = fit_parameters()
  f0 = reference_gaussian.at_d_star_sq(0)
  results = []
  previous_d_min = params.start_interval.d_max + 1
  existing_gaussian = xray_scattering.gaussian([],[])
  while (existing_gaussian.n_ab() < params.max_n_terms):
    if (previous_d_min == params.start_interval.d_min):
      print "%s: Full fit with %d terms. Search stopped." % (
        label, existing_gaussian.n_ab())
      print
      break
    n_terms = existing_gaussian.n_ab() + 1
    fit = find_d_min_multi(
      reference_gaussian=reference_gaussian,
      existing_gaussian=existing_gaussian,
      target_powers=params.target_powers,
      start_interval=params.start_interval,
      n_trial_points_factor=params.n_trial_points_factor,
      n_start_fractions=params.n_start_fractions,
      max_max_error=params.max_max_error,
      min_interval=params.find_d_min_min_interval)
    if (fit.d_min is None):
      print "Warning: No fit: %s n_terms=%d" % (label, n_terms)
      print
      break
    if (previous_d_min < fit.d_min):
      print "Warning: previous d_min was smaller."
    previous_d_min = fit.d_min
    show_fit_summary(
      "Best fit", label, fit.min.final_diff_gaussian, fit.d_min,
      fit.max_error)
    sampling_points = \
      international_tables_sampling_points_and_value_sigmas_up_to(
        d_min=fit.d_min)
    show_literature_fits(
      label=label,
      n_terms=n_terms,
      reference_gaussian=reference_gaussian,
      sampling_points=sampling_points,
      d_min=fit.d_min,
      e_other=fit.max_error)
    fit.min.final_diff_gaussian.show()
    existing_gaussian = fit.min.final_diff_gaussian
    print
    sys.stdout.flush()
    if (plots_dir):
      write_plots(
        plots_dir=plots_dir,
        label=label+"_%d"%n_terms,
        gaussians=[reference_gaussian, fit.min.final_diff_gaussian],
        d_star_sq=sampling_points.d_star_sq)
    g = fit.min.final_diff_gaussian
    results.append(xray_scattering.fitted_gaussian(
      d_min=fit.d_min, a=g.a(),b=g.b()))
  return results

def run(args=[], params=None, verbose=0):
  timer = user_plus_sys_time()
  if (params is None): params = fit_parameters()
  chunk_n = 1
  chunk_i = 0
  if (len(args) > 0 and len(args[0].split(",")) == 2):
    chunk_n, chunk_i = [int(i) for i in args[0].split(",")]
    args = args[1:]
  plots_dir = "plots"
  if (not os.path.isdir(plots_dir)):
    print "***************************************************************"
    print "No plots because target directory does not exist (mkdir %s)." % \
      plots_dir
    print "***************************************************************"
    print
    plots_dir = None
  if (chunk_n > 1):
    assert plots_dir is not None
  results = {}
  results["fit_parameters"] = params
  i_chunk = 0
  for wk in xray_scattering.wk1995_iterator():
    flag = i_chunk % chunk_n == chunk_i
    i_chunk += 1
    if (not flag):
      continue
    if (len(args) > 0 and wk.label() not in args): continue
    results[wk.label()] = incremental_fits(
      label=wk.label(),
      reference_gaussian=wk.fetch(),
      params=params,
      plots_dir=plots_dir,
      verbose=verbose)
    sys.stdout.flush()
    easy_pickle.dump("fits_%02d.pickle" % chunk_i, results)
  print "CPU time: %.2f seconds" % timer.elapsed()

if (__name__ == "__main__"):
  from cctbx.eltbx.gaussian_fit import run
  run(sys.argv[1:])
