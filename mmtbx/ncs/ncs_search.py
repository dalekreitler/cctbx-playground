from __future__ import division
from scitbx.linalg import eigensystem
from scitbx.array_family import flex
from scitbx.math import superpose
from libtbx.utils import Sorry
from scitbx import matrix
import math
import sys
import StringIO
from mmtbx.refinement.flip_peptide_side_chain import should_be_flipped, \
    flippable_sidechains
from time import time

__author__ = 'Youval, massively rewritten by Oleg'

class NCS_groups_container(object):

  def __init__(self):
    """
    A Container for ncs groups
    Note that the first copy is the master ncs

    Attributes:
      iselections (list of flex.size_t):selection array for the complete ASU
      residue_index_list (list): list of list of matching residues indices
      copies (list of lists):List of lists of the ncs copies chain IDs
      transforms (list of transform objects):
    """
    self.iselections = []
    self.residue_index_list = []
    self.copies = []
    self.transforms = []

class Transform(object):

  def __init__(self,
               rotation = None,
               translation = None,
               serial_num = None,
               coordinates_present = None,
               ncs_group_id = None,
               rmsd = 0):
    """
    Basic transformation properties

    Args:
      rotation : Rotation matrix object
      translation: Translation matrix object
      serial_num : (int) Transform serial number
      coordinates_present: equals 1 when coordinates are presents in PDB file
      ncs_group_id : (int) The group ID of the group containing this transform
      rmsd (float): RMS distance between ncs copies
    """
    self.r = rotation
    self.t = translation
    self.serial_num = serial_num
    self.coordinates_present = bool(coordinates_present)
    self.ncs_group_id = ncs_group_id
    self.rmsd = rmsd

class Score_record(object):

  def __init__(self,score=-10000,origin=(0,0)):
    """
    score object used when aligning sequences

    Attributes:
      score (int)
      consecutive_matches (int): num of consecutive matches in current segment
      match_count (int): number of matching residues
      origin (tuple): (row,col) of the matrix cell we from the previous
        alignment step. Used to trace back optimal alignment
      no_altloc (list of bool): False when residue has alternate location
    """
    self.score = score
    self.consecutive_matches = 0
    self.match_count = 0
    self.gap_penalty = 1
    self.origin = origin

class Chains_info(object):
  """ Container for hierarchy analysis """
  def __init__(self):
    self.res_names = []
    self.resid = []
    self.atom_names = []
    self.atom_selection = []
    self.chains_atom_number = 0
    self.no_altloc = []
    self.center_of_coordinates = None

  def __str__(self):
    from StringIO import StringIO
    assert 0
    res = StringIO()
    print >> res, "res_names:", self.res_names
    print >> res, "self.resid", self.resid
    print >> res, "self.atom_names", self.atom_names
    print >> res, "self.atom_selection", self.atom_selection
    print >> res, "self.chains_atom_number", self.chains_atom_number
    print >> res, "self.no_altloc", self.no_altloc
    print >> res, "self.center_of_coordinates", self.center_of_coordinates
    return res.getvalue()


def find_ncs_in_hierarchy(ph,
                          chains_info=None,
                          chain_max_rmsd=5.0,
                          log=None,
                          chain_similarity_threshold=0.85,
                          residue_match_radius=4.0):
  """
  Find NCS relation in hierarchy

  Args:
    ph (object): hierarchy
    use_minimal_master_ncs (bool): use maximal or minimal common chains
        in master ncs groups
    chain_max_rmsd (float): limit of rms difference chains when aligned together
    residue_match_radius (float): max allow distance difference between pairs of matching
      atoms of two residues
    chain_similarity_threshold (float): min similarity between matching chains

  Return:
    groups_list (list of NCS_groups_container objects)
    group_dict (dict):
      keys: tuple of master chain IDs
      values: NCS_groups_container objects
  """
  if not log: log = sys.stdout
  if chains_info is None:
    chains_info = get_chains_info(ph)
  # Get the list of matching chains
  chain_match_list = search_ncs_relations(
    chains_info=chains_info,
    chain_similarity_threshold=chain_similarity_threshold,
    log=None)
  #
  match_dict = clean_chain_matching(
    chain_match_list=chain_match_list,
    ph=ph,
    chain_max_rmsd=chain_max_rmsd,
    residue_match_radius=residue_match_radius)
  # print "match_dict.keys()", match_dict.keys()
  # print "match_dict"
  # for k, v in match_dict.iteritems():
  #   print "  ", k, list(v[0]), list(v[1])
  # assert 0
  #


  # new, the basic way of processing, by Oleg.
  return ncs_grouping_and_group_dict(match_dict, ph)


def get_rmsds2(master_xyz, copy_xyz, cur_ttg):
  """
  This function is for debugging purposes and not called.
  """
  xyz = cur_ttg[2][0].elems * master_xyz + cur_ttg[2][1]
  # rmsd1 = 0
  # if copy_xyz.size() == xyz.size():
  rmsd1 = copy_xyz.rms_difference(xyz)
  xyz = cur_ttg[2][0].elems * master_xyz + cur_ttg[2][1]
  # rmsd2 = 0
  # if copy_xyz.size() == xyz.size():
  rmsd2 = copy_xyz.rms_difference(xyz)
  # print "rmsds:", rmsd1, rmsd2
  return rmsd1, rmsd2

def get_rmsds(hierarchy, cache, cur_ttg, master, copy):
  """
  This function is for debugging purposes and not called.
  Similar check will be performed later in execution and in case of
  wrong grouping will raise Sorry: bad phil records.
  """
  str_sel_m = "chain "+" or chain ".join(cur_ttg[0]+[master])
  str_sel_c = "chain "+" or chain ".join(cur_ttg[1]+[copy])
  sel1 = cache.selection("chain "+" or chain ".join(cur_ttg[0]+[master]))
  sel2 = cache.selection("chain "+" or chain ".join(cur_ttg[1]+[copy]))
  # print "sel1, sel2", str_sel_m, "|", str_sel_c
  master_xyz = hierarchy.select(sel1).atoms().extract_xyz()
  copy_xyz = hierarchy.select(sel2).atoms().extract_xyz()
  xyz = cur_ttg[2][0].elems * master_xyz + cur_ttg[2][1]
  rmsd1 = 0
  if copy_xyz.size() == xyz.size():
    rmsd1 = copy_xyz.rms_difference(xyz)

  str_sel_m = "chain "+" or chain ".join(cur_ttg[0]+[copy])
  str_sel_c = "chain "+" or chain ".join(cur_ttg[1]+[master])
  # print "sel1, sel2", str_sel_m, "|", str_sel_c
  sel1 = cache.selection("chain "+" or chain ".join(cur_ttg[0]+[copy]))
  sel2 = cache.selection("chain "+" or chain ".join(cur_ttg[1]+[master]))
  # print "sel1, sel2", sel1, sel2
  master_xyz = hierarchy.select(sel1).atoms().extract_xyz()
  copy_xyz = hierarchy.select(sel2).atoms().extract_xyz()
  xyz = cur_ttg[2][0].elems * master_xyz + cur_ttg[2][1]
  rmsd2 = 0
  if copy_xyz.size() == xyz.size():
    rmsd2 = copy_xyz.rms_difference(xyz)
  return rmsd1, rmsd2


def get_bool_selection_to_keep(big_selection, small_selection):
  """
  given 2 iselections (they are sorted), returns bool selection of size
  big selection showing what are the matches with small selection.
  Rather fast algorithm but may be beneficial to transfer to C++
  O(n+m), where n,m - sizes of selections
  """
  assert big_selection.size >= small_selection.size()
  result = flex.bool(big_selection.size(), False)
  i_in_big = 0
  i_in_small = 0
  size_small = small_selection.size()
  size_big = big_selection.size()
  n_matches = 0
  nw = 0
  while (i_in_big < size_big) and (i_in_small < size_small):
    if big_selection[i_in_big] == small_selection[i_in_small]:
      result[i_in_big] = True
      i_in_big += 1
      i_in_small += 1
      n_matches += 1
    elif big_selection[i_in_big] > small_selection[i_in_small]:
      i_in_small += 1
      nw += 1
    else:
      i_in_big += 1
  # this assert is optional, in general case it is not guaranteed that
  # all numbers from small selection are present in big selection.
  assert n_matches == size_small, "%d %d" % (n_matches, size_small)
  return result

def get_preliminary_ncs_groups(match_dict):
  pairs = sorted(match_dict.keys())
  chains_in_groups = []
  preliminary_ncs_groups = []
  while len(pairs) > 0:
    # print "  pairs", pairs
    # take the first one, should be new group
    n_not_in_groups = 0
    n_not_in_groups += pairs[0][0] not in chains_in_groups
    n_not_in_groups += pairs[0][1] not in chains_in_groups
    # print "n_not_in_groups", n_not_in_groups
    if n_not_in_groups == 2:
      # make new group
      preliminary_ncs_groups.append({
          pairs[0][0]:pairs[0],
          pairs[0][1]:pairs[0]})
      chains_in_groups.append(pairs[0][0])
      chains_in_groups.append(pairs[0][1])
      curr_masters = pairs[0]
      pairs.pop(0)
      # print "  curr_masters", curr_masters
      # check all the rest pairs to see if they can add something to this group
      pairs_to_remove = []
      for pair in pairs:
        # print "    checking", pair
        if pair[0] == curr_masters[0]:
          if pair[1] not in curr_masters:
            # add pair[1]
            # print "      adding 0"
            if pair[1] not in chains_in_groups:
              preliminary_ncs_groups[-1][pair[1]] = pair
              chains_in_groups.append(pair[1])
            pairs_to_remove.append(pair)

        if pair[1] == curr_masters[0]:
          if pair[0] not in curr_masters:
            # print "      adding 1"
            # add pair[1]
            if pair[0] not in chains_in_groups:
              preliminary_ncs_groups[-1][pair[0]] = pair
              chains_in_groups.append(pair[0])
            pairs_to_remove.append(pair)
      for p in pairs_to_remove:
        pairs.remove(p)

    elif n_not_in_groups == 0:
      # print "    popping the first"
      pairs.pop(0)
    elif n_not_in_groups == 1:
      # should never happen
      # print "    n_not_in_groups==1"
      pairs.pop(0)
      # assert 0
    # print "prel_ncs_gr", preliminary_ncs_groups
  return preliminary_ncs_groups


def ncs_grouping_and_group_dict(match_dict, hierarchy):
  """
  The implementation of simplest way to do NCS grouping. Maximum one chain
  in selection.
  Do the job of minimal_master_ncs_grouping/minimal_ncs_operators_grouping
  and build_group_dict.
  """
  group_dict = {}
  preliminary_ncs_groups = get_preliminary_ncs_groups(match_dict)

  # now we need to just transform preliminary_ncs_groups using match_dict
  # into group_dict. This means that for every dict in preliminary_ncs_groups
  # we need to determine master, and find out rot and transl functions for all
  # the rest chains (selections). Master is going to be the first in
  # alphabetical order.

  group_id = 0
  tr_sn = 1
  for prel_gr_dict in preliminary_ncs_groups:
    # print "==============="
    sorted_gr_chains = sorted(prel_gr_dict.keys())

    # master should be the chain with minimal number of selected atoms
    # just to make it easier filter out the rest of chains
    # print "sorted_gr_chains", sorted_gr_chains
    # print "prel_gr_dict", prel_gr_dict
    min_n_atoms = 1e100
    master = None
    for ch in sorted_gr_chains:
      sel, _,_ = get_info_from_match_dict(match_dict, prel_gr_dict[ch], ch)
      if sel.size() < min_n_atoms:
        min_n_atoms = sel.size()
        master = ch
    assert master is not None
    # print "selected master first:", master

    # second option to master selection:
    # let's try to select common chain to be a master. I'm not sure that this
    # will be always possible though
    # also, we should try to determine the smallest selection for the master
    # chain straight away
    all_pairs = prel_gr_dict.values()
    left = set(all_pairs[0])
    # print "left", left
    # print "all_pairs", all_pairs
    for i in all_pairs[1:]:
      left = left & set(i)
    # should be 1 (a lot of chains) or 2 (if there only 2 chains)
    # if len
    if len(left) == 0:
      # means that all something like
      # all_pairs = [('chain C', 'chain E'), ('chain A', 'chain E'),
      #              ('chain A', 'chain C')]
      # any should work then?...

      # master = all_pairs[0][0]
      master = sorted_gr_chains[0]

    # assert len(left) > 0
    # print "left", left
    elif len(left) > 1:
      master = sorted(left)[0]
    else:
      master = left.pop()


    # selecting smallest master key - for no reason actually
    key_with_smallest_selection = None
    len_of_smallest_selection = 1e100
    for ch, key in prel_gr_dict.iteritems():
      # print "ch, master, key:", ch, master, key
      if master in key:
        master_sel, master_res, master_rmsd = get_info_from_match_dict(
                match_dict, key, master)
        if master_sel.size() < len_of_smallest_selection:
          len_of_smallest_selection = master_sel.size()
          key_with_smallest_selection = key
    # print "key_with_smallest_selection, len_of_smallest_selection",key_with_smallest_selection, len_of_smallest_selection
    # print "selected master second:", master

    assert master is not None
    assert master in key_with_smallest_selection, "%s, %s" % (master, key_with_smallest_selection)

    #
    # Let's do intersection of all master selection to determine
    # the minimum selection suitable to all copies.
    min_master_selection = None
    for ch, key in prel_gr_dict.iteritems():
      if master in key:
        master_sel, master_res, master_rmsd = get_info_from_match_dict(
                match_dict, key, master)
        if min_master_selection is None:
          min_master_selection = master_sel
        else:
          min_master_selection = min_master_selection.intersection(master_sel)
    # print "size of min_master_selection", min_master_selection.size()

    #
    #
    # create a new group
    new_ncs_group = NCS_groups_container()
    tr = Transform(
        rotation=matrix.sqr([1,0,0,0,1,0,0,0,1]),
        translation=matrix.col([0,0,0]),
        serial_num=tr_sn,
        coordinates_present=True,
        ncs_group_id=group_id,
        rmsd=0)
    tr_sn += 1
    # master_sel, master_res, master_rmsd = get_info_from_match_dict(
    #     match_dict,key_with_smallest_selection, master)
    new_ncs_group.iselections.append([min_master_selection])
    new_ncs_group.residue_index_list.append([master_res])
    new_ncs_group.copies.append([master])
    new_ncs_group.transforms.append(tr)

    for ch_copy in sorted_gr_chains:
      master_size = min_master_selection.size()
      copy_sel, copy_res, m_sel = get_copy_master_selections_from_match_dict(
          match_dict, prel_gr_dict, master, ch_copy)
      if copy_sel is None:
        # print " Continue"
        continue
      new_copy_sel = copy_sel
      new_master_sel = min_master_selection

      if copy_sel.size() > min_master_selection.size():
        # clean copy sel
        # print "copy is bigger", copy_sel.size(), min_master_selection.size()
        # print "sizes:", master_sel.size(), m_sel.size()
        filter_sel = get_bool_selection_to_keep(
            big_selection=m_sel,
            small_selection=min_master_selection)
        new_copy_sel = copy_sel.select(filter_sel)
      elif copy_sel.size() < min_master_selection.size():
        # clean master sel and all other copies...
        # should never be the case anymore
        # print "master is bigger", copy_sel.size(), master_sel.size()
        # print "sizes:", master_sel.size(), m_sel.size()
        # print "master:", list(master_sel)
        assert 0
        # filter_sel = get_bool_selection_to_keep(
        #     big_selection=master_sel,
        #     small_selection=m_sel)
        # # print list(filter_sel)
        # new_master_sel = master_sel.select(filter_sel)
        # # print "len new_master_sel", len(new_master_sel)
        # for i in range(len(new_ncs_group.iselections)):
        #   # print "new_ncs_group.iselections", new_ncs_group.iselections
        #   new_ncs_group.iselections[i] = [new_ncs_group.iselections[i][0].select(filter_sel)]
        # master_sel = new_master_sel
        # master_size = master_sel.size()

      r,t,copy_rmsd = my_get_rot_trans(
          ph=hierarchy,
          master_selection=new_master_sel,
          copy_selection=new_copy_sel)
      tr = Transform(
          rotation=r,
          translation=t,
          serial_num=tr_sn,
          coordinates_present=True,
          ncs_group_id=group_id,
          rmsd=copy_rmsd)
      assert master_size == new_copy_sel.size(), "%d %d" % (master_size, new_copy_sel.size())
      new_ncs_group.iselections.append([new_copy_sel])
      new_ncs_group.residue_index_list.append([copy_res])
      new_ncs_group.copies.append([ch_copy])
      new_ncs_group.transforms.append(tr)
      # print "  appended new copy:", ch_copy
      tr_sn += 1
    group_dict[tuple(master)] = new_ncs_group
    master_size = new_ncs_group.iselections[0][0].size()
    for isel_arr in new_ncs_group.iselections[1:]:
      assert master_size ==isel_arr[0].size(), "%d %d" % (master_size, isel_arr[0].size().size())

    # print "new_ncs_group.ise", new_ncs_group.iselections
    # for isele_arr in new_ncs_group.iselections:
    #   print "final selections are:", list(isele_arr[0])
    # print "new_ncs_group.copies", new_ncs_group.copies
    # print "new_ncs_group.residue_index_list", new_ncs_group.residue_index_list
    group_id += 1

  # print "group_dict", group_dict
  # STOP()
  return group_dict


def get_info_from_match_dict(match_dict, key, chain):
  # print "    chain, key in get_info:", chain, key
  assert chain in key, "Mismatch between key and chain %s %s" % (chain, key)
  [sel_1,sel_2,res_1,res_2,_,_,rmsd] = match_dict[key]
  # print "sel_1,sel_2,res_1,res_2,_,_,rmsd", sel_1,sel_2,res_1,res_2,rmsd
  if chain == key[0]:
    return sel_1, res_1, rmsd
  else:
    return sel_2, res_2, rmsd

def get_copy_master_selections_from_match_dict(
    match_dict, prel_gr_dict, master, ch_copy):
  # copy_sel, copy_res, copy_rmsd = get_info_from_match_dict(
  #     match_dict,prel_gr_dict[ch_copy], ch_copy if ch_copy1 is None else ch_copy1)
  # in prel_gr_dict we want to find value with both master and ch_copy
  # return copy_sel, copy_res, m_sel
  key = None
  for v in prel_gr_dict.itervalues():
    if v == (master, ch_copy) or v == (ch_copy, master):
      key = v
      break
  if key is None:
    # print "  key is None, master, ch_copy", master, ch_copy
    return None, None, None
  # print "  key:", key
  [sel_1,sel_2,res_1,res_2,_,_,rmsd] = match_dict[key]
  if master == key[0]:
    return sel_2, res_2, sel_1
  else:
    return sel_1, res_1, sel_2


def make_flips_if_necessary_torsion(const_h, flip_h):
  """ 3 times faster than other procedure."""

  const_h.reset_atom_i_seqs()
  flip_h.reset_atom_i_seqs()
  assert const_h.atoms().size() == flip_h.atoms().size()
  flipped_other_selection = flex.size_t([])
  ch_const = const_h.only_model().chains()
  ch_flip = flip_h.only_model().chains()

  for ch_c, ch_f in zip(ch_const, ch_flip):
    for residue, res_flip in zip(ch_c.residues(), ch_f.residues()):
      if (residue.resname in flippable_sidechains
          and should_be_flipped(residue, res_flip)):
        fl_atom_list = flippable_sidechains[residue.resname]
        iseqs = [0]*residue.atoms().size()
        for i, a in enumerate(residue.atoms()):
          try:
            ind = fl_atom_list.index(a.name)
            if ind == 3 or ind == 5:
              iseqs[i+1] = a.i_seq
            elif ind == 4 or ind == 6:
              iseqs[i-1] = a.i_seq
            else:
              iseqs[i] = a.i_seq
          except ValueError:
            iseqs[i] = a.i_seq
        for i in iseqs:
          flipped_other_selection.append(i)
      else:
        for a in residue.atoms():
          flipped_other_selection.append(a.i_seq)
  # print "flipped_other_selection", list(flipped_other_selection)
  assert flipped_other_selection.size() == const_h.atoms().size()
  return flipped_other_selection

# def make_flips_if_necessary(const_h, flip_h):
#   """ 3 times slower than make_flips_if_necessary_torsion."""
#   def dist(a,b):
#     return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)
#   assert const_h.atoms().size() == flip_h.atoms().size()
#   # this check takes quater of the runtime.
#   if not const_h.contains_protein():
#     return None
#   asc = const_h.atom_selection_cache()
#   # We do alignment only for main chain atoms so that flipped side chains
#   # do not bias it.
#   sel = asc.selection("name N or name CA or name C or name O")
#   ref_sites = const_h.atoms().extract_xyz()
#   other_sites = flip_h.atoms().extract_xyz()
#   ref_sites_for_fit = const_h.atoms().extract_xyz().select(sel)
#   other_sites_for_fit = flip_h.atoms().extract_xyz().select(sel)
#   lsq_fit_obj = superpose.least_squares_fit(
#     reference_sites = ref_sites_for_fit,
#     other_sites     = other_sites_for_fit)
#   r = lsq_fit_obj.r
#   t = lsq_fit_obj.t
#   flip_sites_best = r.elems*other_sites + t.elems
#   # first reset i_seq to use them directly with sites arrays
#   # they should not contain alt confs
#   # XXX multiple chains when there is HETATM with a ligand after TER.
#   flipped_other_selection = flex.size_t([])
#   const_h.reset_atom_i_seqs()
#   flip_h.reset_atom_i_seqs()
#   for ch in const_h.only_model().chains():
#     for residue in ch.only_conformer().residues():
#       if (residue.resname in ["GLU", "ASP", "PHE", "HIS", "LEU",
#                               "ASN", "GLN", "ARG", "VAL", "TYR"] and
#           residue.atoms().size() > 1):
#         # find interesting pair and decide on flip straight away
#         flippable_iseqs = []
#         atoms = residue.atoms()
#         i = 0
#         while i < len(atoms):
#           iname = atoms[i].name.strip()
#           i1name = atoms[i+1].name.strip() if i < len(atoms)-1 else ""
#           if (len(iname) == len(i1name) == 3 and
#               iname[-2] == i1name[-2] and
#               iname[-1] == "1" and i1name[-1] == "2"):
#             flippable_iseqs.append((atoms[i].i_seq, atoms[i+1].i_seq))
#             is1 = atoms[i].i_seq
#             is2 = atoms[i+1].i_seq
#             dist_non_flip = dist(ref_sites[is1], flip_sites_best[is1]) + \
#                 dist(ref_sites[is2], flip_sites_best[is2])
#             dist_flip = dist(ref_sites[is1], flip_sites_best[is2]) + \
#                 dist(ref_sites[is2], flip_sites_best[is1])
#             # print "dist non flip, flip:", dist_non_flip, dist_flip
#             if dist_flip < dist_non_flip - 0.5:
#               flipped_other_selection.append(atoms[i+1].i_seq)
#               flipped_other_selection.append(atoms[i].i_seq)
#             else:
#               flipped_other_selection.append(atoms[i].i_seq)
#               flipped_other_selection.append(atoms[i+1].i_seq)
#             i += 1
#           else:
#             flipped_other_selection.append(atoms[i].i_seq)
#           i += 1
#
#       else:
#         # residue is not flippable, goes straight to flipped_other_selection
#         for a in residue.atoms():
#           flipped_other_selection.append(a.i_seq)
#   assert flipped_other_selection.size() == const_h.atoms().size()
#   # print "flipped_other_selection", list(flipped_other_selection)
#   return flipped_other_selection


def clean_chain_matching(chain_match_list,ph,
                         chain_max_rmsd=10.0,
                         residue_match_radius=4.0):
  """
  Remove all bad matches from chain_match_list

  Args:
    ph (object): hierarchy
    chain_match_list (list): list of
      [chain_ID_1, chain_ID_2, sel_1, sel_2,res_m/res_c similarity]
      chain_ID (str), sel_1/2 (list of lists)
      res_m/res_c (lists): indices of the aligned components
      similarity (float): similarity between chains
    chain_max_rmsd (float): limit of rms difference chains
    residue_match_radius (float): max allow distance difference between pairs of matching
      atoms of two residues
    chain_similarity_threshold (float): min similarity between matching chains

  Returns:
    match_dict(dict): key:(chains_id_a,chains_id_b)
                      val:[selection_a,selection_b,
                           res_list_a,res_list_b,rot,trans,rmsd]
  """
  # remove all non-matching pairs, where similarity == 0
  match_list = [x for x in chain_match_list if x[4] > 0]
  match_dict = {}
  # print "match_list", match_list
  for match in match_list:
    [ch_a_id,ch_b_id,list_a,list_b,res_list_a,res_list_b,similarity] = match
    t0 = time()
    sel_a = make_selection_from_lists(list_a)
    sel_b = make_selection_from_lists(list_b)

    other_h = ph.select(sel_a)
    other_atoms = other_h.atoms()
    ref_h = ph.select(sel_b)
    ref_atoms = ref_h.atoms()
    #
    # Here we want to flip atom names, even before chain alignment, so
    # we will get correct chain RMSD

    # flipped_other_selection = make_flips_if_necessary(ref_h.deep_copy(), other_h.deep_copy())
    flipped_other_selection = make_flips_if_necessary_torsion(
        ref_h.deep_copy(), other_h.deep_copy())
    # if flipped_other_selection is not None:
    other_sites = other_atoms.select(flipped_other_selection).extract_xyz()
    # else:
    #   other_sites = other_atoms.extract_xyz()
    ref_sites = ref_atoms.extract_xyz()
    lsq_fit_obj = superpose.least_squares_fit(
      reference_sites = ref_sites,
      other_sites     = other_sites)
    r = lsq_fit_obj.r
    t = lsq_fit_obj.t
    # todo: find r_2*A = r*A + t (where the translation is zero)
    # use B = r*A + t, r_2*A = B , r_2 = B*A.inverse()
    other_sites_best = lsq_fit_obj.other_sites_best_fit()
    rmsd = round(ref_sites.rms_difference(other_sites_best),4)
    # print "chain rmsd after flip:", rmsd
    if rmsd <= chain_max_rmsd:
      # get the chains atoms and convert selection to flex bool
      sel_aa,sel_bb,res_list_a,res_list_b,ref_sites,other_sites_best = \
        remove_far_atoms(
          list_a, list_b,
          res_list_a,res_list_b,
          ref_sites,lsq_fit_obj.other_sites_best_fit(),
          residue_match_radius=residue_match_radius)
      if sel_a.size() > 0:
        match_dict[ch_a_id,ch_b_id]=[sel_aa,sel_bb,res_list_a,res_list_b,r,t,rmsd]
  return match_dict

def remove_far_atoms(list_a, list_b,
                     res_list_a,res_list_b,
                     ref_sites,other_sites,
                     residue_match_radius=4.0):
  """
  When comparing lists of matching atoms, remove residues where some atoms are
  are locally misaligned, for example when matching residues are
  perpendicular to each other rather than being close to parallel.

  The criteria used:
  For each matching residues, the difference between distance of farthest
  matching atoms pair and the distance of closest pair mast be < residue_match_radius

  Args:
    list_a, list_a (list of list): list of residues atoms
    res_list_a,res_list_b (list): list of residues in chains
    ref_sites,other_sites (flex.vec3): atoms coordinates
    residue_match_radius (float): max allow distance difference

  Returns:
    Updated arguments:
      sel_a,sel_b,
      res_list_a_new,res_list_b_new,
      ref_sites_new,other_sites_new
  """
  # check every residue for consecutive distance
  # print "list_a"
  # print list(list_a[0])
  # print "list_b", list(list_b)
  # print "res_list_a", res_list_a
  # print "res_list_b", res_list_b
  res_list_a_new = []
  res_list_b_new = []
  ref_sites_new = flex.vec3_double([])
  other_sites_new = flex.vec3_double([])
  sel_a = flex.size_t([])
  sel_b = flex.size_t([])
  current_pos = 0
  for i in xrange(len(res_list_a)):
    # find the matching atoms form each residue (work on small sections)
    res_len = list_a[i].size()
    res_ref_sites = ref_sites[current_pos:current_pos+res_len]
    res_other_sites = other_sites[current_pos:current_pos+res_len]
    current_pos += res_len
    xyz_diff = abs(res_ref_sites.as_double() - res_other_sites.as_double())
    (min_d,max_d,_) = xyz_diff.min_max_mean().as_tuple()
    # print "current match radius:", max_d-min_d
    if (max_d - min_d) <= residue_match_radius:
      ref_sites_new.extend(res_ref_sites)
      other_sites_new.extend(res_other_sites)
      sel_a.extend(list_a[i])
      sel_b.extend(list_b[i])
      res_list_a_new.append(res_list_a[i])
      res_list_b_new.append(res_list_b[i])
    else:
      pass
      # print "removing poorly matching residue:",i,max_d - min_d
  return sel_a,sel_b,res_list_a_new,res_list_b_new,ref_sites_new,other_sites_new


def find_same_transform(r,t,transforms):
  """

  Not used.

  Check if the rotation r and translation t exist in the transform dictionary.
  Note that there can be both inverse and regular match. Return the
  non-transpose if exist.

  Args:
    r (matrix.sqr): rotation
    t (matrix.col): translation
    transforms (dict): dictionary of all transforms

  Returns:
    tr_num: (str) transform serial number
    is_transpose: (bool) True if the matching transform is transpose
  """

  is_transpose = False
  tr_num = None
  for k,v in transforms.iteritems():
    if hasattr(v,'r'):
      rr = v.r
      tt = v.t
    else:
      (rr,tt) = v[2]
    is_the_same, is_transpose_flag = is_same_transform(r, t, rr, tt)
    if is_the_same:
      if is_transpose_flag:
        # when transpose is found, keep it but continue the search
        tr_num = k
        is_transpose  = True
      else:
        # found non-transform match
        return k, False
  return tr_num, is_transpose

def get_sequence_from_array(arr):
  from iotbx.pdb import amino_acid_codes
  aa_3_as_1 = amino_acid_codes.one_letter_given_three_letter
  res = ""
  for r in arr:
    res += aa_3_as_1.get(r)
  return res

def search_ncs_relations(ph=None,
                         chains_info = None,
                         chain_similarity_threshold=0.85,
                         log=None):
  """
  Search for NCS relations between chains or parts of chains, in a protein
  hierarchy

  Args:
    ph (object): hierarchy
    chains_info (dict): values are object containing
      chains (str): chain IDs OR selections string
      res_name (list of str): list of residues names
      resid (list of str): list of residues sequence number, resid
      atom_names (list of list of str): list of atoms in residues
      atom_selection (list of list of list of int): the location of atoms in ph
      chains_atom_number (list of int): list of number of atoms in each chain

  Returns:
    msg (str): message regarding matching residues with different atom number
    chain_match_list (list): list of
      [chain_ID_1,chain_ID_2,sel_1,sel_2,res_sel_m, res_sel_c,similarity]
      chain_ID (str), sel_1/2 (flex.size_t),
      res_sel_m/c (lists of lists): indices of the aligned components
      similarity (float): similarity between chains
    We use sel_2 to avoid problems when residues have different number of atoms
  """
  # print "searching ncs relations..."
  if not log: log = StringIO.StringIO()
  if not chains_info:
    assert bool(ph)
    chains_info = get_chains_info(ph)
  # collect all chain IDs
  chain_match_list = []
  msg = ''
  sorted_ch = sorted(chains_info)

  n_chains = len(sorted_ch)
  chains_in_copies = set()
  for i in xrange(n_chains-1):
    m_ch_id = sorted_ch[i]

    if m_ch_id in chains_in_copies:
      continue

    master_n_res = len(chains_info[m_ch_id].res_names)
    seq_m = chains_info[m_ch_id].res_names
    if master_n_res == 0:
      continue
    # get residue lists for master
    for j in xrange(i+1,n_chains):
      c_ch_id = sorted_ch[j]
      copy_n_res = len(chains_info[c_ch_id].res_names)
      frac_d = min(copy_n_res,master_n_res)/max(copy_n_res,master_n_res)
      if frac_d < chain_similarity_threshold:
        if (chain_similarity_threshold == 1):
          msg = 'NCS copies are not identical'
          break
        else:
          # print "Strange exit"
          continue
      seq_c = chains_info[c_ch_id].res_names
      # get residue lists for copy
      res_sel_m, res_sel_c, similarity = mmtbx_res_alignment(
          seq_a=seq_m,seq_b=seq_c,
          min_percent=chain_similarity_threshold)
      sel_m, sel_c,res_sel_m,res_sel_c,new_msg = get_matching_atoms(
        chains_info,m_ch_id,c_ch_id,res_sel_m,res_sel_c)
      msg += new_msg
      if res_sel_m:
        # add only non empty matches
        rec = [m_ch_id,c_ch_id,sel_m,sel_c,res_sel_m,res_sel_c,similarity]
        chain_match_list.append(rec)
      # Collect only very good matches, to allow better similarity search
      if similarity > chain_similarity_threshold:
        chains_in_copies.add(c_ch_id)
        # print "  good"
  # loop over all chains
  if msg:
    print >> log,msg
  if (chain_similarity_threshold == 1) and msg:
    # must be identical
    raise Sorry('NCS copies are not identical')
  return chain_match_list

def mmtbx_res_alignment(seq_a, seq_b,
                        min_percent=0.85, atomnames=False):
  # Check for the basic cases (shortcut for obvious cases)
  a = len(seq_a)
  b = len(seq_b)
  if (a == 0) or (b == 0): return [], [], 0
  if seq_a == seq_b: return range(a), range(a), 1.0
  norm_seq_a = seq_a
  norm_seq_b = seq_b
  if not atomnames:
    norm_seq_a = ""
    norm_seq_b = ""
    from iotbx.pdb.amino_acid_codes import one_letter_given_three_letter, \
        one_letter_given_three_letter_modified_aa
    merged_one_given_three = one_letter_given_three_letter.copy()
    merged_one_given_three.update(one_letter_given_three_letter_modified_aa)
    merged_one_given_three.update({
        "  A": "A",
        "  C": "C",
        "  G": "G",
        "  U": "U",
        " DA": "A",
        " DC": "C",
        " DG": "G",
        " DT": "T"})
    for l in seq_a:
      one_letter = merged_one_given_three.get(l, 'X')
      norm_seq_a += one_letter
    for l in seq_b:
      one_letter = merged_one_given_three.get(l, 'X')
      norm_seq_b += one_letter
  from mmtbx.alignment import align
  obj = align(
      norm_seq_a,
      norm_seq_b,
      gap_opening_penalty=1, # default
      gap_extension_penalty=0.5, # default is 1
      similarity_function="identity")
  alignment = obj.extract_alignment()
  sim1 = alignment.calculate_sequence_identity()
  # print "Sequence identity is", sim1
  # alignment.pretty_print(block_size=60)
  al_a, al_b = alignment.exact_match_selections()
  # alignment.pretty_print()

  if sim1 < min_percent:
    # chains are to different, return empty arrays
    return flex.size_t([]), flex.size_t([]), 0
  return al_a, al_b, sim1


def get_matching_atoms(chains_info,a_id,b_id,res_num_a,res_num_b):
  """
  Get selection of matching chains, match residues atoms
  We keep only residues with continuous matching atoms

  Residues with alternative locations and of different size are excluded

  Args:
    chains_info (object): object containing
      chains (str): chain IDs or selections string
      res_name (list of str): list of residues names
      resid (list of str): list of residues sequence number, resid
      atom_names (list of list of str): list of atoms in residues
      atom_selection (list of list of list of int): the location of atoms in ph
      chains_atom_number (list of int): list of number of atoms in each chain
    a_id,b_id (str): Chain IDs
    res_num_a/b (list of int): indices of matching residues position

  Returns:
    sel_a/b (list of lists): matching atoms selection
    res_num_a/b (list of int): updated res_num_a/b
    msg (str): message regarding matching residues with different atom number
  """
  sel_a = []
  sel_b = []
  # check if any of the residues has alternate locations
  a_altloc = bool(chains_info[a_id].no_altloc)
  if a_altloc:
    a_altloc = chains_info[a_id].no_altloc.count(False) > 0
  b_altloc = bool(chains_info[b_id].no_altloc)
  if b_altloc:
    b_altloc = chains_info[b_id].no_altloc.count(False) > 0
  test_altloc = a_altloc or b_altloc
  #
  res_num_a_updated = []
  res_num_b_updated = []
  residues_with_different_n_atoms = []
  for (i,j) in zip(res_num_a,res_num_b):
    # iterate over atoms in residues
    sa = flex.size_t(chains_info[a_id].atom_selection[i])
    sb = flex.size_t(chains_info[b_id].atom_selection[j])
    dif_res_size = sa.size() != sb.size()
    atoms_names_a = chains_info[a_id].atom_names[i]
    atoms_names_b = chains_info[b_id].atom_names[j]
    resid_a = chains_info[a_id].resid[i]
    force_check_atom_order = dif_res_size
    altloc = False
    if test_altloc:
      if a_altloc:
        altloc |= (not chains_info[a_id].no_altloc[i])
      if b_altloc:
        altloc |= (not chains_info[b_id].no_altloc[j])
    if force_check_atom_order:
      # select only atoms that exist in both residues
      atoms_a,atoms_b,similarity = mmtbx_res_alignment(
        seq_a=atoms_names_a, seq_b=atoms_names_b,
        min_percent=0.2, atomnames=True)
      # get the number of the atom in the chain
      sa = flex.size_t(atoms_a) + sa[0]
      sb = flex.size_t(atoms_b) + sb[0]
    if dif_res_size or altloc:
      residues_with_different_n_atoms.append(resid_a)
      if altloc:
        sa = flex.size_t([])
        sb = flex.size_t([])
    # keep only residues with continuous matching atoms
    if sa.size() != 0:
      res_num_a_updated.append(i)
      res_num_b_updated.append(j)
      sel_a.append(sa)
      sel_b.append(sb)
  if residues_with_different_n_atoms:
    problem_res_nums = {x.strip() for x in residues_with_different_n_atoms}
    msg = "NCS related residues with different number of atoms, selection"
    msg += a_id + ':' + b_id + '\n['
    msg += ','.join(problem_res_nums) + ']\n'
  else:
    msg = ''
  return sel_a,sel_b,res_num_a_updated,res_num_b_updated,msg

def make_selection_from_lists(sel_list):
  """ Convert a list of lists to flex.size_t selection array  """
  sel_list_extended = [x for y in sel_list for x in y]
  sel_set = set(sel_list_extended)
  assert len(sel_list_extended) == len(sel_set)
  sel_list_extended.sort()
  return flex.size_t(sel_list_extended)

def get_chains_info(ph, selection_list=None):
  """
  Collect information about chains or segments of the hierarchy according to
  selection strings
  Exclude water atoms
  When there are alternate conformations, we use the first one

  Args:
    ph : pdb_hierarchy

  Returns:
    chains_info (dict): values are object containing
      chains (str): chain IDs OR selections string
      res_name (list of str): list of residues names
      resid (list of str): list of residues sequence number, resid
      atom_names (list of list of str): list of atoms in residues
      atom_selection (list of list of list of int): the location of atoms in ph
      chains_atom_number (list of int): list of number of atoms in each chain
    exclude_water (bool): exclude water
  """

  chains_info =  {}
  # asc = ph.atom_selection_cache()
  model  = ph.models()[0]
  # build chains_info from hierarchy
  for ch in model.chains():
    if not chains_info.has_key(ch.id):
      chains_info[ch.id] = Chains_info()
      # This is very time-consuming
      # ph_sel = ph.select(asc.selection("chain '%s'" % ch.id))
      # coc = flex.vec3_double([ph_sel.atoms().extract_xyz().mean()])
      # chains_info[ch.id].center_of_coordinates = coc
      chains_info[ch.id].center_of_coordinates = None
    chains_info[ch.id].chains_atom_number += ch.atoms().size()
    resids = chains_info[ch.id].resid
    res_names = chains_info[ch.id].res_names
    atom_names = chains_info[ch.id].atom_names
    atom_selection = chains_info[ch.id].atom_selection
    no_altloc = chains_info[ch.id].no_altloc
    conf = ch.conformers()[0]
    len_conf = len(ch.conformers())
    for res in conf.residues():
      x = res.resname
      resids.append(res.resid())
      res_names.append(x)
      atoms = res.atoms()
      atom_names.append(list(atoms.extract_name()))
      atom_selection.append(list(atoms.extract_i_seq()))
      no_altloc.append(res.is_pure_main_conf or len_conf==1)
    chains_info[ch.id].resid = resids
    chains_info[ch.id].res_names = res_names
    chains_info[ch.id].atom_names = atom_names
    chains_info[ch.id].atom_selection = atom_selection
    chains_info[ch.id].no_altloc = no_altloc
  return chains_info


def inverse_transform(r,t):
  """ inverse rotation and translation """
  r = r.transpose()
  t = - r*t
  return r,t

def angle_between_rotations(v1,v2):
  """ get angle between two vectors"""
  cos_angle = v1.dot(v2)
  result = math.acos(min(1,cos_angle))
  result *= 180/math.pi
  return result

def get_rotation_vec(r):
  """ get the eigen vector associated with the eigen value 1"""
  eigen = eigensystem.real_symmetric(r.as_sym_mat3())
  eigenvectors = eigen.vectors()
  eigenvalues = eigen.values()
  i = list(eigenvalues.round(4)).index(1)
  return eigenvectors[i:(i+3)]

def is_same_transform(r1,t1,r2,t2):
  """
  Check if transform is the same by comparing rotations and the result of
  applying rotation and translation on
  a test vector

  Args:
    r1, r2: Rotation matrices
    t1, t2: Translation vectors

  Returns:
    (bool,bool) (is_the_same, is_transpose)
  """
  # Allowed deviation for values and angle
  eps=0.1
  angle_eps=5.0
  if (not r1.is_zero()) and (not r2.is_zero()):
    assert r1.is_r3_rotation_matrix(rms_tolerance=0.001)
    assert r2.is_r3_rotation_matrix(rms_tolerance=0.001)
    # test vector
    xyz = flex.vec3_double([(11,103,523),(-500.0,2.0,10.0),(0.0,523.0,-103.0)])
    a_ref = (r1.elems * xyz + t1).as_double()
    rt, tt = inverse_transform(r1,t1)
    a_ref_transpose = (rt.elems * xyz + tt).as_double()
    v1 = get_rotation_vec(r1)
    v2 = get_rotation_vec(r2)
    a = (r2.elems * xyz + t2).as_double()
    d = (a_ref-a)
    d = (d.dot(d))**.5/a.size()
    dt = (a_ref_transpose-a)
    dt = (dt.dot(dt))**.5/a.size()
    ang = angle_between_rotations(v1,v2)
    d_ang = min(ang, (180 - ang))
    if (d_ang < angle_eps) and (d < eps):
      return True, False
    elif (d_ang < angle_eps) and (dt < eps):
      return True, True
    else:
      return False, False
  else:
    return False, False


def my_get_rot_trans(
    ph,
    master_selection,
    copy_selection):
  """
  Get rotation and translation using superpose.

  This function is used only when phil parameters are provided. In this case
  we require the selection of NCS master and copies to be correct.
  Correct means:
    1) residue sequence in master and copies is exactly the same
    2) the number of atoms in master and copies is exactly the same

  One can get exact selection strings by ncs_object.show(verbose=True)

  Args:
    ph : hierarchy
    master/copy_selection: master and copy iselections
  """

  atoms = ph.atoms()
  # master
  other_sites = atoms.select(master_selection).extract_xyz()
  # copy
  ref_sites = atoms.select(copy_selection).extract_xyz()
  assert other_sites.size() == ref_sites.size(), "%d, %d" % (
      other_sites.size(), ref_sites.size())
  if ref_sites.size() > 0:
    lsq_fit_obj = superpose.least_squares_fit(
        reference_sites = ref_sites,
        other_sites     = other_sites)
    r = lsq_fit_obj.r
    t = lsq_fit_obj.t
    rmsd = ref_sites.rms_difference(lsq_fit_obj.other_sites_best_fit())
    return r,t,rmsd
  else:
    assert 0, "strange thing..."


def get_rot_trans(ph,
                  master_selection,
                  copy_selection,
                  chain_max_rmsd=0.02):
  """
  Get rotation and translation using superpose.

  This function is used only when phil parameters are provided. In this case
  we require the selection of NCS master and copies to be correct.
  Correct means:
    1) residue sequence in master and copies is exactly the same
    2) the number of atoms in master and copies is exactly the same

  One can get exact selection strings by ncs_object.show(verbose=True)

  Args:
    ph : pdb.hierarchy
    master/copy_selection (str): master and copy selection strings
    chain_max_rmsd (float): limit of rms difference between chains to be considered
      as copies

  Returns:
    r: rotation matrix
    t: translation vector
    rmsd (float): RMSD between master and copy
    msg (str): error messages
  """
  msg = ''
  r_zero = matrix.sqr([0]*9)
  t_zero = matrix.col([0,0,0])
  #
  if ph:
    cache = ph.atom_selection_cache().selection
    master_ncs_ph = ph.select(cache(master_selection))
    ncs_copy_ph = ph.select(cache(copy_selection))
    seq_m,res_ids_m  = get_residue_sequence(master_ncs_ph)
    seq_c,res_ids_c = get_residue_sequence(ncs_copy_ph)
    res_sel_m, res_sel_c, similarity = mmtbx_res_alignment(
        seq_m, seq_c, min_percent=0)
    # res_sel_m, res_sel_c, similarity = res_alignment(
    #   seq_a=seq_m,seq_b=seq_c,
    #   min_contig_length=0,min_percent=0)
    m_atoms = master_ncs_ph.atoms()
    c_atoms = ncs_copy_ph.atoms()
    # Check that master and copy are identical
    if (similarity != 1) or (m_atoms.size() != c_atoms.size()) :
      return r_zero,t_zero,0,'Master and Copy selection do not exactly match'
    # master
    other_sites = m_atoms.extract_xyz()
    # copy
    ref_sites = c_atoms.extract_xyz()
    if ref_sites.size() > 0:
      lsq_fit_obj = superpose.least_squares_fit(
          reference_sites = ref_sites,
          other_sites     = other_sites)
      r = lsq_fit_obj.r
      t = lsq_fit_obj.t
      rmsd = ref_sites.rms_difference(lsq_fit_obj.other_sites_best_fit())
      if rmsd > chain_max_rmsd:
        return r_zero,t_zero,0,msg
    else:
      return r_zero,t_zero,0,'No sites to compare.\n'
    return r,t,round(rmsd,4),msg
  else:
    return r_zero,t_zero,0,msg


def get_residue_sequence(ph):
  """
  Get a list of residues numbers and names from hierarchy "ph", excluding
  water molecules

  Args:
    ph (hierarchy): hierarchy of a single chain

  Returns:
    res_list_new (list of str): list of residues names
    resid_list_new (list of str): list of residues number
  """
  res_list_new = []
  resid_list_new = []
  for res_info in ph.atom_groups():
    x = res_info.resname
    if x.lower() != 'hoh':
      # Exclude water
      res_list_new.append(x)
      resid_list_new.append(res_info.parent().resseq)
  return res_list_new,resid_list_new
