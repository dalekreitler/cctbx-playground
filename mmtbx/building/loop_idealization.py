from __future__ import division
import iotbx.pdb
import mmtbx.utils
from mmtbx.rotamer.rotamer_eval import RotamerEval
from mmtbx.monomer_library import idealized_aa
from libtbx.utils import Sorry, null_out
from mmtbx.validation import ramalyze
from mmtbx.building.loop_closure.ccd import ccd_python
from mmtbx.building.loop_closure import utils, starting_conformations
from mmtbx.pdbtools import truncate_to_poly_gly
from mmtbx.secondary_structure.build import side_chain_placement, \
    set_xyz_smart
from mmtbx.refinement.geometry_minimization import minimize_wrapper_for_ramachandran

import boost.python
ext = boost.python.import_ext("mmtbx_validation_ramachandran_ext")
from mmtbx_validation_ramachandran_ext import rama_eval

loop_idealization_master_phil_str = """
loop_idealization
{
  enabled = True
    .type = bool
  change_non_rama_outliers = True
    .type = bool
    .help = Allow changing non-outlier ramachandran residues angles
  output_prefix = rama_fixed
    .type = str
  minimize_whole = True
    .type = bool
  force_rama_fixes = True
    .type = bool
    .help = If true, the procedure will pick and apply the best variant even \
      if all of them are above thresholds to be picked straight away. \
      Alternatively, when False, the procedure will accept failure and leave \
      a ramachandran outlier intact.
  save_states = False
    .type = bool
    .help = Save states of CCD. Generates a states file for every model.
  number_of_ccd_trials = 5
    .type = int
    .help = How many times we are trying to fix outliers in the same chain
  variant_search_level = 2
    .type = int
    .help = how thoroughly variants will be explored (1-3)
  variant_number_cutoff = 50
    .type = int
    .help = how many first variants to take from generated
}
"""

master_phil = iotbx.phil.parse(loop_idealization_master_phil_str)

class loop_idealization():
  def __init__(self,
               pdb_hierarchy,
               params=None,
               secondary_structure_annotation=None,
               log=null_out(),
               verbose=True):
    if len(pdb_hierarchy.models()) > 1:
      raise Sorry("Multi-model files are not supported")
    self.original_pdb_h = pdb_hierarchy
    self.secondary_structure_annotation=secondary_structure_annotation
    xrs = pdb_hierarchy.extract_xray_structure()
    asc = pdb_hierarchy.atom_selection_cache()
    self.resulting_pdb_h = pdb_hierarchy.deep_copy()
    self.params = self.process_params(params)
    self.log = log
    self.verbose = verbose
    self.r = rama_eval()
    self.rotamer_manager = RotamerEval()
    ram = ramalyze.ramalyze(pdb_hierarchy=pdb_hierarchy)
    self.p_initial_rama_outliers = ram.out_percent
    self.p_before_minimization_rama_outliers = None
    self.p_after_minimiaztion_rama_outliers = None

    berkeley_count = utils.list_rama_outliers_h(self.resulting_pdb_h).count("\n")
    self.berkeley_p_before_minimization_rama_outliers = \
        berkeley_count/float(self.resulting_pdb_h.overall_counts().n_residues)*100

    # self.berkeley_p_before_minimization_rama_outliers = None
    self.berkeley_p_after_minimiaztion_rama_outliers = None
    self.ref_exclusion_selection = ""
    number_of_ccd_trials = 0
    # print "logic expr outcome:", (number_of_ccd_trials < 10 and self.berkeley_p_before_minimization_rama_outliers > 0.001)
    # print number_of_ccd_trials < 10
    # print self.berkeley_p_before_minimization_rama_outliers > 0.001
    while (number_of_ccd_trials < self.params.number_of_ccd_trials
        and self.berkeley_p_before_minimization_rama_outliers > 0.001):
      print "CCD try number, outliers:", number_of_ccd_trials, self.berkeley_p_before_minimization_rama_outliers
      number_of_ccd_trials += 1
      for chain in self.resulting_pdb_h.only_model().chains():
        print >> self.log, "Idealizing chain %s" % chain.id
        selection = "protein and chain %s and (name N or name CA or name C or name O)" % chain.id
        sel = asc.selection("chain %s" % chain.id)
        chain_h = self.resulting_pdb_h.select(sel)
        m = chain_h.only_model()
        i = 0
        cutted_chain_h = None
        for c in m.chains():
          if i == 0:
            cutted_chain_h = iotbx.pdb.hierarchy.new_hierarchy_from_chain(c)
          else:
            print >> self.log, "WARNING!!! Duplicating chain ids! Only the first chain will be processed."
            print >> self.log, "  Removing chain %s with %d residues" % (c.id, len(c.residues()))
            m.remove_chain(c)
          i += 1
        exclusions, ch_h = self.idealize_chain(
            hierarchy=(cutted_chain_h if cutted_chain_h else chain_h))
        if ch_h is not None:
          print "Setting new coordinates"
          set_xyz_smart(self.resulting_pdb_h, ch_h)
          for resnum in exclusions:
            selection += " and not resseq %s" % resnum
        self.ref_exclusion_selection += "(%s) or " % selection
      #
      # dumping and reloading hierarchy to do proper rounding of coordinates
      self.resulting_pdb_h = iotbx.pdb.input(
          source_info=None,
          lines=self.resulting_pdb_h.as_pdb_string()).construct_hierarchy()
      berkeley_count = utils.list_rama_outliers_h(self.resulting_pdb_h).count("\n")
      self.berkeley_p_before_minimization_rama_outliers = \
          berkeley_count/float(self.resulting_pdb_h.overall_counts().n_residues)*100
    if len(self.ref_exclusion_selection) > 0:
      self.ref_exclusion_selection = self.ref_exclusion_selection[:-3]
    self.resulting_pdb_h.write_pdb_file(file_name="%s_before_minimization.pdb" % self.params.output_prefix)
    ram = ramalyze.ramalyze(pdb_hierarchy=self.resulting_pdb_h)
    self.p_before_minimization_rama_outliers = ram.out_percent

    duke_count = ram.get_outliers_count_and_fraction()[0]
    if berkeley_count != duke_count:
      print >> self.log, "Discrepancy between berkeley and duke after ccd:", berkeley_count, duke_count

    if self.params.minimize_whole:
      print >> self.log, "minimizing whole thing..."
      print >> self.log, "self.ref_exclusion_selection", self.ref_exclusion_selection
      print >> sel
      minimize_wrapper_for_ramachandran(
          hierarchy=self.resulting_pdb_h,
          xrs=xrs,
          original_pdb_h=self.original_pdb_h,
          excl_string_selection=self.ref_exclusion_selection,
          log=None,
          ss_annotation=self.secondary_structure_annotation)
      # self.resulting_pdb_h.write_pdb_file(file_name="%s_all_minized.pdb" % self.params.output_prefix)
      ram = ramalyze.ramalyze(pdb_hierarchy=self.resulting_pdb_h)
      self.p_after_minimiaztion_rama_outliers = ram.out_percent
      berkeley_count = utils.list_rama_outliers_h(self.resulting_pdb_h).count("\n")
      duke_count = ram.get_outliers_count_and_fraction()[0]
      self.berkeley_p_after_minimiaztion_rama_outliers = \
          berkeley_count/float(self.resulting_pdb_h.overall_counts().n_residues)*100
      if berkeley_count != duke_count:
        print >> self.log, "Discrepancy between berkeley and duke after min:", berkeley_count, duke_count
    # return new_h

  def process_params(self, params):
    if params is None:
      params = master_phil.fetch().extract()
      params.loop_idealization.enabled = True
    if hasattr(params, "model_idealization"):
      p_pars = params.model_idealization
    else:
      assert hasattr(params, "enabled") and hasattr(params, "change_non_rama_outliers"), \
          "Something wrong with parameters passed to model_idealization"
      p_pars = params

    assert isinstance(p_pars.enabled, bool)
    assert isinstance(p_pars.change_non_rama_outliers, bool)
    return p_pars

  def idealize_chain(self, hierarchy):
    # check no ac:
    for c in hierarchy.chains():
      if len(c.conformers()) > 1:
        raise Sorry("Alternative conformations are not supported.")
      if "UNK" in c.get_residue_names_padded():
        raise Sorry("UNK residues are not supported.")
    working_h = hierarchy.deep_copy()
    working_h.reset_atom_i_seqs()
    rama_results = []
    ranges_for_idealization = []
    print >> self.log, "rama outliers for input hierarchy:"
    rama_out_resnums = self.get_resnums_of_chain_rama_outliers(
        working_h)
    if len(rama_out_resnums) == 0:
      return None, None
    # get list of residue numbers that should be excluded from reference
    list_of_reference_exclusion = []
    for resnum in rama_out_resnums:
      start_rn, prev_rn = get_res_nums_around(hierarchy, resnum, 2, 2)
      list_of_reference_exclusion += [start_rn, resnum, prev_rn]
    out_i = 0
    for rama_out_resnum in rama_out_resnums:
      print >> self.log
      print >> self.log, "Fixing outlier:", rama_out_resnum
      self.log.flush()
      new_h = self.fix_rama_outlier(
        pdb_hierarchy=working_h,
        out_res_num=rama_out_resnum,
        prefix=self.params.output_prefix,
        minimize=False)
      print >> self.log, "listing outliers after loop minimization"
      outp = utils.list_rama_outliers_h(new_h, self.r)
      print >> self.log, outp
      self.log.flush()
      fn = "%s_after_loop_%d.pdb" % (self.params.output_prefix, out_i)
      # print >> self.log, "  writing file %s" % fn
      # new_h.write_pdb_file(file_name=fn)
      working_h = new_h
      out_i += 1
    return list_of_reference_exclusion, new_h

  def ccd_solution_is_ok(self,
      anchor_rmsd, mc_rmsd, ccd_radius, change_all_angles, change_radius):
    adaptive_mc_rmsd = {1:3.0, 2:3.5, 3:4.0}
    if (mc_rmsd < adaptive_mc_rmsd[ccd_radius] and anchor_rmsd < 0.3):
      return True
    elif ccd_radius == 3 and change_all_angles and change_radius == 2:
      # we are desperate and trying the most extensive search,
      # this deserves relaxed criteria...
      return mc_rmsd < 5 and anchor_rmsd < 0.4



  def fix_rama_outlier(self,
      pdb_hierarchy, out_res_num, prefix="", minimize=True):

    original_pdb_h = pdb_hierarchy.deep_copy()
    chain_id = original_pdb_h.only_model().only_chain().id
    all_results = []
    variants_searches = [
        ((1, False, 0),1),
        ((2, False, 0),1),
        ((3, False, 0),2),
        ((2, True, 1),1),
        ((3, True, 1),2),
        ((3, True, 2),3),
    ]
    decided_variants = []
    for variant, level in variants_searches:
      if level <= self.params.variant_search_level:
        decided_variants.append(variant)

    for ccd_radius, change_all, change_radius in decided_variants:
    # while ccd_radius <= 3:
      print >> self.log, "  Starting optimization with radius, change_all, change_radius:", ccd_radius, change_all, change_radius
      self.log.flush()
      moving_h, moving_ref_atoms_iseqs, fixed_ref_atoms = get_fixed_moving_parts(
          pdb_hierarchy=pdb_hierarchy,
          out_res_num=out_res_num,
          n_following=ccd_radius,
          n_previous=ccd_radius)
      moving_h_set = None
      if change_all:
        moving_h_set = starting_conformations.get_all_starting_conformations(
            moving_h,
            change_radius,
            cutoff=self.params.variant_number_cutoff,
            # log=self.log,
            )
      else:
        moving_h_set = starting_conformations.get_starting_conformations(
            moving_h,
            cutoff=self.params.variant_number_cutoff,
            # log=self.log,
            )

      if len(moving_h_set) == 0:
        # outlier was fixed before somehow...
        # or there's a bug in get_starting_conformatiosn
        return original_pdb_h

      for i, h in enumerate(moving_h_set):
        ccd_obj = ccd_python(fixed_ref_atoms, h, moving_ref_atoms_iseqs)
        ccd_obj.run()
        resulting_rmsd = ccd_obj.resulting_rmsd
        states = ccd_obj.states
        n_iter = ccd_obj.n_iter
        if self.params.save_states:
          states.write(file_name="%s%s_%d_%s_%d_%i_states.pdb" % (chain_id, out_res_num, ccd_radius, change_all, change_radius, i))

        # resulting_rmsd, states, n_iter = ccd(
        #     fixed_ref_atoms, h, moving_ref_atoms_iseqs, moving_h)

        mc_rmsd = get_main_chain_rmsd_range(moving_h, h, all_atoms=True)
        print >> self.log, "Resulting anchor and backbone RMSDs, n_iter for model %d:" % i,
        print >> self.log, resulting_rmsd, ",", mc_rmsd, ",", n_iter
        self.log.flush()
        #
        # setting new coordinates
        #
        moved_with_side_chains_h = pdb_hierarchy.deep_copy()
        set_xyz_smart(moved_with_side_chains_h, h)
        #
        # placing side-chains
        #
        # moved_with_side_chains_h.write_pdb_file(
        #     file_name="%s_before_sc_placement_%d.pdb" % (prefix, i))
        placing_range = get_res_nums_around(moved_with_side_chains_h,
            center_resnum=out_res_num,
            n_following=ccd_radius,
            n_previous=ccd_radius,
            include_intermediate=True)
        place_side_chains(moved_with_side_chains_h, original_pdb_h,
            self.rotamer_manager, placing_range)
        # moved_with_side_chains_h.write_pdb_file(
        #     file_name="%s_after_sc_placement_%d.pdb" % (prefix, i))


        #
        # finalizing with geometry_minimization
        #
        # !!! This is the condition of acceptance of transformation!
        # if mc_rmsd < adaptive_mc_rmsd[ccd_radius]:
        all_results.append((moved_with_side_chains_h.deep_copy(), mc_rmsd, resulting_rmsd, n_iter))
        if self.ccd_solution_is_ok(
            anchor_rmsd=resulting_rmsd,
            mc_rmsd=mc_rmsd,
            ccd_radius=ccd_radius,
            change_all_angles=change_all,
            change_radius=change_radius):
          print "Choosen result (mc_rmsd, anchor_rmsd, n_iter):", mc_rmsd, resulting_rmsd, n_iter
          self.log.flush()
          if minimize:
            print >> self.log, "minimizing..."
            moved_with_side_chains_h.write_pdb_file(
                file_name="%s_result_before_min_%d.pdb" % (prefix, i))
            minimize_wrapper_for_ramachandran(
                hierarchy=moved_with_side_chains_h,
                xrs=xrs,
                original_pdb_h=original_pdb_h,
                log=self.log,
                ss_annotation=self.secondary_structure_annotation)
          moved_with_side_chains_h.write_pdb_file(
              file_name="%s_result_minimized_%d.pdb" % (prefix, i))
          final_rmsd = get_main_chain_rmsd_range(moved_with_side_chains_h,
              original_pdb_h, placing_range)
          print >> self.log, "FINAL RMSD after minimization:", final_rmsd
          return moved_with_side_chains_h


    all_results.sort(key=lambda tup: tup[1])
    print "ALL RESULTS:"
    i = 0
    for ar in all_results:
      print ar[1:],
      if ar[2] < 0.4:
        fn = "variant_%d.pdb" % i
        ar[0].write_pdb_file(file_name=fn)
        print fn
        i += 1
      else:
        print "  no output"
    if self.params.force_rama_fixes:
      # find and apply the best varian from all_results. This would be the one
      # with the smallest rmsd given satisfactory closure
      print >> self.log, "Applying the best found variant:",
      i = 0
      while i < len(all_results) and all_results[i][2] > 0.1:
        i += 1
      # apply
      # === duplication!!!!
      if i < len(all_results):
        print >> self.log, all_results[i][1:]
        if minimize:
          print >> self.log, "minimizing..."
          all_results[i][0].write_pdb_file(
              file_name="%s_result_before_min_%d.pdb" % (prefix, i))
          minimize_wrapper_for_ramachandran(
              hierarchy=all_results[i][0],
              xrs=xrs,
              original_pdb_h=original_pdb_h,
              log=self.log,
              ss_annotation=self.secondary_structure_annotation)
        all_results[i][0].write_pdb_file(
            file_name="%s_result_minimized_%d.pdb" % (prefix, i))
        final_rmsd = get_main_chain_rmsd_range(all_results[i][0],
            original_pdb_h, placing_range)
        print >> self.log, "FINAL RMSD after minimization:", final_rmsd
        return moved_with_side_chains_h
      else:
        print >> self.log, " NOT FOUND!"
        for i in all_results:
          print >> self.log, i[1:]
      # === end of duplication!!!!


    else:
      print >> self.log, "Epic FAIL: failed to fix rama outlier"
      print >> self.log, "  Options were: (mc_rmsd, resultign_rmsd, n_iter)"
      for i in all_results:
        print >> self.log, i[1:]
    # STOP()
    return original_pdb_h

  def get_resnums_of_chain_rama_outliers(self, pdb_hierarchy):
    phi_psi_atoms = utils.get_phi_psi_atoms(pdb_hierarchy)
    # print "len phi psi atoms", len(phi_psi_atoms)
    result = []
    rama_results = []
    ranges_for_idealization = []
    # print >> self.log, "rama outliers for input hierarchy:"
    list_of_reference_exclusion = []
    outp = utils.list_rama_outliers_h(pdb_hierarchy, self.r)
    print >> self.log, outp
    for phi_psi_pair, rama_key in phi_psi_atoms:
      # print "resseq:", phi_psi_pair[0][2].parent().parent().resseq
      ev = utils.rama_evaluate(phi_psi_pair, self.r, rama_key)
      # print "  ev", ev
      rama_results.append(ev)
      if ev == ramalyze.RAMALYZE_OUTLIER:
        resnum = phi_psi_pair[0][2].parent().parent().resseq
        result.append(resnum)
    return result

def place_side_chains(hierarchy, original_h,
    rotamer_manager, placing_range):
  ideal_res_dict = idealized_aa.residue_dict()
  asc = original_h.atom_selection_cache()
  gly_atom_names = set([" N  ", " CA ", " C  ", " O  "])
  for rg in hierarchy.residue_groups():
    if rg.resseq in placing_range:
      # cut extra atoms
      ag = rg.only_atom_group()
      for atom in ag.atoms():
        if (atom.name not in gly_atom_names):
          ag.remove_atom(atom=atom)
      # get ag from original hierarchy
      orig_ag = original_h.select(asc.selection("resseq %d" % rg.resseq_as_int())
          ).models()[0].chains()[0].residue_groups()[0].atom_groups()[0]
      # get ideal
      ideal_ag = ideal_res_dict[ag.resname.lower()].models()[0].chains()[0].\
        residue_groups()[0].atom_groups()[0]
      # print "got to placement"
      side_chain_placement(ag, orig_ag, rotamer_manager)

def get_res_nums_around(pdb_hierarchy, center_resnum, n_following, n_previous,
    include_intermediate=False):
  residue_list = list(
      pdb_hierarchy.only_model().only_chain().only_conformer().residues())
  center_index = None
  for i in range(len(residue_list)):
    if residue_list[i].resseq == center_resnum:
      center_index = i
      break
  # print "start/end resids", residue_list[i-n_previous].resseq, residue_list[i+n_following].resseq
  # print "center i, len", center_index, len(residue_list)
  if not include_intermediate:
    return residue_list[max(0,center_index-n_previous)].resseq, \
        residue_list[min(len(residue_list)-1,center_index+n_following)].resseq
  else:
    res = []
    for i in range(max(0,center_index-n_previous),
        min(len(residue_list)-1,center_index+n_following+1)):
      res.append(residue_list[i].resseq)
    return res


def get_fixed_moving_parts(pdb_hierarchy, out_res_num, n_following, n_previous):
  # limitation: only one  chain in pdb_hierarchy!!!
  original_pdb_h = pdb_hierarchy.deep_copy()
  start_res_num, end_res_num = get_res_nums_around(
      pdb_hierarchy, out_res_num, n_following, n_previous)

  xrs = original_pdb_h.extract_xray_structure()
  truncate_to_poly_gly(pdb_hierarchy, start_res_num, end_res_num)
  cache = pdb_hierarchy.atom_selection_cache()
  # print "selectioin:", "resid %d through %d" % (start_res_num, end_res_num)
  m_selection = cache.selection("(name N or name CA or name C or name O) and resid %s through %s" % (start_res_num, end_res_num))
  moving_h = pdb_hierarchy.select(m_selection)
  moving_h.reset_atom_i_seqs()
  # print dir(moving_h)
  # STOP()
  m_cache = moving_h.atom_selection_cache()
  # print "len inp h atoms", pdb_hierarchy.atoms().size()
  # print "len moving_h atoms", moving_h.atoms().size()
  moving_ref_atoms_iseqs = []
  # here we need N, CA, C atoms from the end_res_num residue
  eff_end_resnum = end_res_num
  sel = m_cache.selection("resid %s" % end_res_num)
  while len(moving_h.select(sel).atoms()) == 0:
    eff_end_resnum -= 1
    sel = m_cache.selection("resid %s" % eff_end_resnum)

  sel = m_cache.selection("resid %s and name N" % eff_end_resnum)
  a = moving_h.select(sel).atoms()[0]
  moving_ref_atoms_iseqs.append(a.i_seq)
  fixed_N = a.detached_copy()

  sel = m_cache.selection("resid %s and name CA" % eff_end_resnum)
  a = moving_h.select(sel).atoms()[0]
  moving_ref_atoms_iseqs.append(a.i_seq)
  fixed_CA = a.detached_copy()

  sel = m_cache.selection("resid %s and name C" % eff_end_resnum)
  a = moving_h.select(sel).atoms()[0]
  moving_ref_atoms_iseqs.append(a.i_seq)
  fixed_C = a.detached_copy()

  fixed_ref_atoms = [fixed_N, fixed_CA, fixed_C]

  return moving_h, moving_ref_atoms_iseqs, fixed_ref_atoms

def get_main_chain_rmsd_range(
    hierarchy, original_h, all_atoms=False, placing_range=None):
  rmsd = 0
  mc_atoms = None
  if all_atoms:
    mc_atoms = ["N", "CA", "C", "O"]
  else:
    mc_atoms = ["N", "CA", "C"]
  for m_atom, ref_atom in zip(hierarchy.atoms(), original_h.atoms()):
    if m_atom.name.strip() in mc_atoms:
      if (placing_range is None or
          m_atom.parent().parent().resseq in placing_range):
        rmsd += m_atom.distance(ref_atom)**2
  return rmsd**0.5
