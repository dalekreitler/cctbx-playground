import cctbx.xray.structure_factors
from cctbx.xray import ext
from cctbx.xray.structure import structure as cctbx_xray_structure
from cctbx.array_family import flex
import scitbx.lbfgs
from scitbx.python_utils.misc import adopt_init_args

class lbfgs:

  def __init__(self, target_functor, gradient_flags, xray_structure,
                     lbfgs_termination_params=None,
                     lbfgs_core_params=None,
                     cos_sin_table=0001,
                     structure_factor_algorithm=None):
    adopt_init_args(self, locals())
    self.structure_factors_from_scatterers = \
      cctbx.xray.structure_factors.from_scatterers(
        miller_set=self.target_functor.f_obs(),
        cos_sin_table=cos_sin_table)
    self.structure_factor_gradients = \
      cctbx.xray.structure_factors.gradients(
        miller_set=self.target_functor.f_obs(),
        cos_sin_table=cos_sin_table)
    self.n = xray_structure.n_parameters(gradient_flags)
    self.x = flex.double(self.n, 0)
    self._scatterers_start = xray_structure.scatterers()
    self._scattering_dict = xray_structure.scattering_dict()
    self._d_min = self.target_functor.f_obs().d_min()
    self.first_target_value = None
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      core_params=lbfgs_core_params)
    self.apply_shifts()
    del self._scatterers_start
    del self._scattering_dict
    del self._d_min
    self.compute_target(compute_gradients=00000)
    self.final_target_value = self.target_result.target()

  def apply_shifts(self):
    scatterers_shifted = ext.minimization_apply_shifts(
      self.xray_structure.unit_cell(),
      self.xray_structure.space_group_info().type(),
      self._scatterers_start,
      self._scattering_dict,
      self.gradient_flags,
      self.x,
      self._d_min)
    self.xray_structure.replace_scatterers(scatterers_shifted)

  def compute_target(self, compute_gradients):
    self.f_calc = self.structure_factors_from_scatterers(
      xray_structure=self.xray_structure,
      miller_set=self.target_functor.f_obs(),
      algorithm=self.structure_factor_algorithm).f_calc()
    self.target_result = self.target_functor(
      self.f_calc,
      compute_gradients)

  def __call__(self):
    if (self.first_target_value is None):
      assert self.x.all_eq(0)
    else:
      self.apply_shifts()
    self.compute_target(compute_gradients=0001)
    self.f = self.target_result.target()
    if (self.first_target_value is None):
      self.first_target_value = self.f
    sf = self.structure_factor_gradients(
      xray_structure=self.xray_structure,
      miller_set=self.target_functor.f_obs(),
      d_target_d_f_calc=self.target_result.derivatives(),
      gradient_flags=self.gradient_flags,
      n_parameters=self.x.size(),
      algorithm=self.structure_factor_algorithm)
    self.g = sf.packed()
    return self.x, self.f, self.g
