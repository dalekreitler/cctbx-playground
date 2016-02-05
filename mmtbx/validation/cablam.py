from __future__ import division
# (jEdit options) :folding=explicit:collapseFolds=1:
from mmtbx.validation import residue, validation, atom
from cctbx import geometry_restraints
from iotbx import pdb, file_reader
from mmtbx.rotamer.n_dim_table import NDimTable #handles contours
from libtbx import easy_pickle #NDimTables are stored as pickle files
import libtbx.load_env
import libtbx.phil.command_line
import os, sys

#{{{ global constants
#-------------------------------------------------------------------------------
MAINCHAIN_ATOMS = [" N  "," C"  ," O  "," CA "]
CA_PSEUDOBOND_DISTANCE = 4.5
PEPTIDE_BOND_DISTANCE = 2.0
#2.0A is the peptide bond cutoff used by O for model building and potentially-
#  fragmented chains. O's generous cutoff seemed appropriate since I expect to
#  process in-progress models with this program
#RLab (probably) uses 1.4 +- 0.3, official textbook is about 1.33
#O's Ca-Ca distance is 4.5A
#Tom T suggests 2.5A as low-end cutoff for Ca-Ca distance
#--He also suggests using the pdb interpreter's internal chain break def

CABLAM_OUTLIER_CUTOFF = 0.01
CABLAM_DISFAVORED_CUTOFF = 0.05
CA_GEOM_CUTOFF = 0.005
#These cutoffs are used to identify outliers in CaBLAM parameter spaces
#These values were set heuristically through inspection of known outliers

ALPHA_CUTOFF = 0.001
BETA_CUTOFF = 0.0001
THREETEN_CUTOFF = 0.001
#These cutoffs are used to identify probable secondary structure residues
#Individual secondary structure residues are assembled into complete secondary
#  structure elements
#These values were set heuristically through inspection of known secondary
#  structure in low-resolution models
#-------------------------------------------------------------------------------
#}}}

#{{{ phil
#-------------------------------------------------------------------------------
#commandline parameters
master_phil = libtbx.phil.parse("""
cablam {
  pdb_infile = None
    .type = path
    .help = '''input PDB file'''
  output = *text kin full_kin records records_and_pdb oneline
    .type = choice
    .help = '''choose output type:
    =text for default colon-separated residue-by-residue validation
    =kin for outlier markup in kinemage format
    =full_kin for outlier markup appended to structure - opens in KiNG
    =records for PDB-style HELIX/SHEET records
    =records_and_pdb for PDB-style HELIX/SHEET records attached to a PDB file
    =oneline for a one-line structure summary
    '''
  outliers_only = False
    .type = bool
    .help = '''compress certain outputs to show only outliers'''
  help = False
    .type = bool
    .help = '''help and data interpretation messages'''
}
""",process_includes=True)
#-------------------------------------------------------------------------------
#}}}

#{{{ usage
#-------------------------------------------------------------------------------
def usage():
  #prints help text
  sys.stderr.write("""
phenix.cablam file.pdb [options ...]

Options:

  pdb_infile=filename   path to input PDB file
                          structure coordinates file readable by Phenix.
                          Supports .pdb, .ent, .cif, etc.

  output=(choose one)   text : default output.  Prints machine-readable
                          columnated and colon-separated validation text to
                          screen.
                        kin : prints kinemage markup for validation to screen
                        full_kin : prints kinemage markup and struture kinamge
                          to screen
                        records : prints pdb-style HELIX and SHEET records to
                          screen, based on CaBLAM's identification of secondary
                          structure
                        records_and_pdb : prints pdb-style HELIX and SHEET
                          records to screen, followed by PDB file formatted
                          coordinates for the input structure
                        oneline : prints single-line summary of CaBLAM
                          validation statistics

  outliers_only=False   compresses certain outputs (text) to show only outlier
                          residues

  help=False            prints this usage text, plus an interpretation guide
""")
#-------------------------------------------------------------------------------
#}}}

#{{{ interpretation
#-------------------------------------------------------------------------------
def interpretation():
  #prints a brief guide to interpreting CaBLAM results
  sys.stderr.write("""
---------------- *** CaBLAM validation data interpretation *** -----------------

CaBLAM uses CA-trace geometry to validate protein backbone. Since the CA trace
  is more reliably modeled from low-resolution density than the full backbone
  trace, CaBLAM can identify modeling errors and intended secondary structure
  elements in low-resolution models.

Text:
  Text output is provided in colon-separated format. The columns are as follows:
    residue : A residue identifier
    outlier_type : If the residue is an outlier by CaBLAM metrics, the type will
      appear here. There are 3 types of outliers:
        CaBLAM Outlier is similar to Ramachandran 'Outlier' in severity
        CaBLAM Disfavored is similar to Ramachandran 'Allowed'
        CA Geom Outlier indicates a severe error in the CA trace
    contour_level : The CaBLAM contour percentile score for the residue.
      0.05 (5%) or lower is Disfavored.  0.01 (1%) or lower is Outlier.
    ca_contour_level : The CA geometry contour percentile score for the residue.
      Serves as a sanity check for CaBLAM's other metrics.
      0.005 (0.5%) or lower is Outlier.
    sec struc recommendation : Secondary structure identification for this
      residue. There are 3 possible identifications - try alpha helix,
      try beta sheet, and try three-ten.
    alpha score : The alpha helix contour percentile score for this residue.
      0.001 (0.1%) or better is probable helix.
    beta score : The beta strand contour percentile score for this residue.
      0.0001 (0.01%) or better is probable strand.
    three-ten score : The 3-10 helix contour percentile score for this residue.
      0.001 (0.1%) or better is probable helix.

Kinemage:
  Kinemage output is available for visual markup of structures in KiNG.
""")
#-------------------------------------------------------------------------------
#}}}

#{{{ math functions
#-------------------------------------------------------------------------------
def perptersect(a1, a2, b1):
  #Finds the line from a1 to a2, drops a perpendicular to it from b1, and returns
  #  the point of intersection.
  A = [a2[0]-a1[0], a2[1]-a1[1], a2[2]-a1[2]]
  #Find the slope of line A in each direction, A is in vector notation
  t = (A[0]*(b1[0]-a1[0]) + A[1]*(b1[1]-a1[1]) + A[2]*(b1[2]-a1[2])) / ((A[0]**2)+(A[1]**2)+(A[2]**2))
  #Solve the parametric equations (dot of perpendiculars=0). . .
  b2 = [a1[0]+A[0]*t, a1[1]+A[1]*t, a1[2]+A[2]*t]
  # . . . and use the result to find the new point b2 on the line
  return b2

def calculate_mu(CA1,CA2,CA3,CA4):
  #dihedral calculation for CA trace
  if None in [CA1,CA2,CA3,CA4]:
    return None
  return geometry_restraints.dihedral(sites=[CA1.xyz,CA2.xyz,CA3.xyz,CA4.xyz],
    angle_ideal=180, weight=1).angle_model

def calculate_ca_virtual_angle(CA1,CA2,CA3):
  #angle calculation for CA trace
  if None in [CA1,CA2,CA3]:
    return None
  return geometry_restraints.angle(sites=[CA1.xyz,CA2.xyz,CA3.xyz],
    angle_ideal=120, weight=1).angle_model

def calculate_nu(CA1,CA2,CA3,O1,O2):
  #dihedral calculation for peptide plane orientations
  if None in [CA1,CA2,CA3,O1,O2]:
    return None
  X1 = perptersect(CA1.xyz,CA2.xyz,O1.xyz)
  X2 = perptersect(CA2.xyz,CA3.xyz,O2.xyz)
  return geometry_restraints.dihedral(sites=[O1.xyz,X1,X2,O2.xyz],
    angle_ideal=180, weight=1).angle_model

def calculate_omega(CA1,C1,N2,CA2):
  #dihedral calculation for peptide bond
  if None in [CA1,C1,N2,CA2]:
    return None
  return geometry_restraints.dihedral(sites=[CA1.xyz,C1.xyz,N2.xyz,CA2.xyz],
    angle_ideal=180, weight=1).angle_model
#-------------------------------------------------------------------------------
#}}}

#{{{ contour fetching
#-------------------------------------------------------------------------------
def fetch_peptide_expectations():
  #This function finds, unpickles, and returns N-Dim Tables of expected residue
  #  behavior for use in determining cablam outliers
  #The return object is a dict keyed by residue type:
  #  'general','gly','transpro','cispro'
  #This set of contours defines peptide behavior (mu_in,mu_out,nu)
  categories = ['general','gly','transpro','cispro']
  unpickled = {}
  for category in categories:
    picklefile = libtbx.env.find_in_repositories(
      relative_path=(
        "chem_data/cablam_data/cablam.8000.expected."+category+".pickle"),
      test=os.path.isfile)
    if (picklefile is None):
      sys.stderr.write("\nCould not find a needed pickle file for category "+
        category+" in chem_data.\nExiting.\n")
      sys.exit()
    ndt = easy_pickle.load(file_name=picklefile)
    unpickled[category] = ndt
  return unpickled

def fetch_ca_expectations():
  #This function finds, unpickles, and returns N-Dim Tables of expected residue
  #  behavior for use in determining ca geometry outliers
  #The return object is a dict keyed by residue type:
  #  'general','gly','transpro','cispro'
  #This set of contours defines CA trace quality (mu_in,mu_d_out,ca_virtual)
  categories = ['general','gly','transpro','cispro']
  unpickled = {}
  for category in categories:
    picklefile = libtbx.env.find_in_repositories(
      relative_path=(
        "chem_data/cablam_data/cablam.8000.expected."+category+"_CA.pickle"),
      test=os.path.isfile)
    if (picklefile is None):
      sys.stderr.write("\nCould not find a needed pickle file for category "+
        category+" in chem_data.\nExiting.\n")
      sys.exit()
    ndt = easy_pickle.load(file_name=picklefile)
    unpickled[category] = ndt
  return unpickled

def fetch_motif_contours():
  #This function finds, unpickles, and returns N-Dim Tables of secondary structure
  #  structure behavior for use in determining what an outlier residue might
  #  really be
  #The return object is a dict keyed by secondary structure type:
  #  'regular_beta','loose_alpha','loose_threeten'
  motifs = ['regular_beta','loose_alpha','loose_threeten']
  unpickled = {}
  for motif in motifs:
    picklefile = libtbx.env.find_in_repositories(
      relative_path=(
        "chem_data/cablam_data/cablam.8000.motif."+motif+".pickle"),
      test=os.path.isfile)
    if (picklefile is None):
      sys.stderr.write("\nCould not find a needed pickle file for motif "+
        motif+" in chem_data.\nExiting.\n")
      sys.exit()
    ndt = easy_pickle.load(file_name=picklefile)
    unpickled[motif] = ndt
  return unpickled
#-------------------------------------------------------------------------------
#}}}

#{{{ cablam data storage classes
#-------------------------------------------------------------------------------
class cablam_geometry():
  #holds cablam's geometry parameters for one residue
  def __init__(self, mu_in=None, mu_out=None, nu=None, ca_virtual=None, omega=None):
    self.mu_in = mu_in
    self.mu_out = mu_out
    self.nu = nu
    self.ca_virtual = ca_virtual
    self.omega=omega

class cablam_score():
  #holds cablam contour scores for one residue
  def __init__(self,cablam=None,c_alpha_geom=None,alpha=None,beta=None,threeten=None):
    self.cablam = cablam
    self.c_alpha_geom = c_alpha_geom
    self.alpha = alpha
    self.beta = beta
    self.threeten = threeten

class cablam_feedback():
  #holds outliers status and secondary structure identifications for one residue
  def __init__(self):
    self.cablam_outlier = None
    self.cablam_disfavored = None
    self.c_alpha_geom_outlier = None
    self.alpha=None
    self.beta=None
    self.threeten=None

#cablam results are stored by chains and by conformers within chains, in
#  parallel with the conformer organization of the pdb hierarchy.  These classes
#  handle that organization of the results
class cablam_chain():
  #chain-level organization for cablam results
  def __init__(self):
    self.conf_names = []
    self.confs = {}

class cablam_conf():
  #conformer-level organization for cablam results
  def __init__(self):
    self.conf_name = None
    self.results = {}
    self.sec_struc_records = []

class secondary_structure_segment():
  #holds a secondary structure element identified by cablam
  def __init__(self, start, end, segment_type, segment_length):
    self.start = start
    self.end = end
    self.segment_type = segment_type
    self.segment_length = segment_length
#-------------------------------------------------------------------------------
#}}}

#{{{ cablam_result class
#-------------------------------------------------------------------------------
class cablam_result(residue):
  """
  Result class for cablam analysis
  """
  __slots__ = residue.__slots__ + [
    "residue",
    "prevres",
    "nextres",
    "alt",
    "chain",
    "_outlier_i_seqs",
    "has_ca",
    "has_mc",
    "has_any_alts",
    "has_mc_alts",
    "mc_alts",
    #"is_outlier",
    "alts",
    "measures",
    "scores",
    "feedback"
    ]

  #{{{ mp_id
  #-----------------------------------------------------------------------------
  def mp_id(self):
    #Returns an id consistent with MolProbity 'cnit' ids
    #Formatted as: ccnnnnilttt
    #  c: 2-char Chain ID, space for none
    #  n: sequence number, right justified, space padded
    #  i: insertion code, space for none
    #  l: alternate ID, space for none
    #  t: residue type (ALA, LYS, etc.), all caps left justified, space padded
    ##id_str = self.residue.id_str()
    ##| A  75 |
    ##chain = id_str[0:2]
    chain = self.chain
    resnum = self.residue.resseq
    ins = self.residue.icode
    resname = self.residue.resname
    alt = self.alt
    if alt == '':
      alt = ' '
    return chain + resnum + ins + alt + resname
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ sorting_id
  #-----------------------------------------------------------------------------
  def sorting_id(self):
    #Returns an id used for sorting residues
    #Formatted as: ccnnnni
    #  c: 2-char Chain ID, space for none
    #  n: sequence number, right justified, space padded
    #  i: insertion code, space for none
    chain = self.chain
    resnum = self.residue.resseq
    ins = self.residue.icode
    return chain + resnum + ins
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ get_atom
  #-----------------------------------------------------------------------------
  #returns the desired atom from a hierarchy residue object
  def get_atom(self, atom_name):
    for atom in self.residue.atoms():
      if atom.name == atom_name: return atom
    else: return None
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ get_cablam_atoms
  #-----------------------------------------------------------------------------
  def get_cablam_atoms(self):
    #finds and returns all the atoms necessary for the CaBLAM calculations for
    #  this residue
    #returned atoms are:
    #res0_CA
    #res1_CA, res1_O, res1_C
    #res2_CA, res2_O, res2_N (res2 is the current residue)
    #res3_CA
    #res4_CA
    #returns None for missing/unfound atoms
    atom_set = {}
    atom_set["res2_CA"]=self.get_atom(" CA ")
    atom_set["res2_O"]= self.get_atom(" O  ")#for nu
    atom_set["res2_N"]= self.get_atom(" N  ")#for omega
    if self.prevres:
      atom_set["res1_CA"]=self.prevres.get_atom(" CA ")
      atom_set["res1_O"]= self.prevres.get_atom(" O  ")#for nu
      atom_set["res1_C"]= self.prevres.get_atom(" C  ")#for omega
      if self.prevres.prevres:
        atom_set["res0_CA"]=self.prevres.prevres.get_atom(" CA ")
      else:
        atom_set["res0_CA"]=None
    else:
      atom_set["res1_CA"]=None
      atom_set["res1_O"]= None
      atom_set["res1_C"]= None
      atom_set["res0_CA"]=None
    if self.nextres:
      atom_set["res3_CA"]=self.nextres.get_atom(" CA ")
      if self.nextres.nextres:
        atom_set["res4_CA"]=self.nextres.nextres.get_atom(" CA ")
      else:
        atom_set["res4_CA"]=None
    else:
      atom_set["res3_CA"]=None
      atom_set["res4_CA"]=None
    return atom_set
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ check_atoms
  #-----------------------------------------------------------------------------
  #checks whether atoms necessary for cablam calculations are present in this
  #  residue.
  #Sets has_ca = True if the residue has a CA atom (minimum cablam requirement)
  #Sets has_mc = True if the residue has all 4 mainchain heavy atoms
  def check_atoms(self):
    if self.get_atom(' CA ') is None: pass
    else:
      self.has_ca = True
      for atom_name in [' N  ',' C  ',' O  ']:
        if self.get_atom(atom_name) is None:
          break
      else:
        self.has_mc = True
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ link_residues
  #-----------------------------------------------------------------------------
  def link_residues(self, previous_result):
    #CaBLAM calculations depend on traversing sequential residues both backwards
    #  and forwards in sequence.  This function creates "links" between
    #  sequential residues.
    #Links are established if the residues are within probable bonding distance.
    #prevres and nextres values default to None during cablam_data.__init__()
    if previous_result is None:
      #no previous residue to link to
      return #Default link is None
    elif not previous_result.has_ca or not self.has_ca:
      #previous residue's proximity cannot be checked
      return
    elif not previous_result.has_mc or not self.has_mc:
      #CA-trace-only: use CA to check proximity
      ca1 = previous_result.get_atom(' CA ').xyz
      ca2 = self.get_atom(' CA ').xyz
      cadist = ((ca1[0]-ca2[0])**2 + (ca1[1]-ca2[1])**2 + (ca1[2]-ca2[2])**2)**0.5
      if cadist > CA_PSEUDOBOND_DISTANCE:
        #CA atoms are too far apart
        return
      else:
        #
        previous_result.nextres = self
        self.prevres = previous_result
    else: #has full set of mc heavy atoms available
      # use previous ' C  ' and current ' N  ' to check proximity
      c = previous_result.get_atom(' C  ').xyz
      n = self.get_atom(' N  ').xyz
      peptidedist = ((c[0]-n[0])**2 + (c[1]-n[1])**2 + (c[2]-n[2])**2)**0.5
      if peptidedist > PEPTIDE_BOND_DISTANCE:
        #atoms that would peptide bond are too far apart
        return
      else:
        previous_result.nextres = self
        self.prevres = previous_result
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ calculate_cablam_geomtery
  #-----------------------------------------------------------------------------
  def calculate_cablam_geometry(self):
    #populates self.measures with geometric measures relevant to CaBLAM
    #these measures are: mu_in, mu_out, nu, ca_virtual, omega
    atom_set = self.get_cablam_atoms()
    self.measures = cablam_geometry(
      mu_in  = calculate_mu(atom_set['res0_CA'],atom_set['res1_CA'],atom_set['res2_CA'],atom_set['res3_CA']),
      mu_out = calculate_mu(atom_set['res1_CA'],atom_set['res2_CA'],atom_set['res3_CA'],atom_set['res4_CA']),
      nu = calculate_nu(atom_set['res1_CA'],atom_set['res2_CA'],atom_set['res3_CA'],atom_set['res1_O'],atom_set['res2_O']),
      ca_virtual = calculate_ca_virtual_angle(atom_set['res1_CA'],atom_set['res2_CA'],atom_set['res3_CA']),
      omega = calculate_omega(atom_set['res1_CA'],atom_set['res1_C'],atom_set['res2_N'],atom_set['res2_CA'])
      )
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ contour_category
  #-----------------------------------------------------------------------------
  def contour_category(self):
    #determines the category of the current residue so that it can be paired
    #  with the correct contours
    #these categories are: 'general', 'gly', 'transpro', 'cispro'
    resname = self.residue.resname.upper()
    if resname == "GLY": return "gly"
    elif resname == "PRO":
      if self.measures.omega < 90 and self.measures.omega > -90:
        return "cispro"
      else: return "transpro"
    else: return "general"
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ calculate_contour_values
  #-----------------------------------------------------------------------------
  def calculate_contour_values(self, cablam_contours, ca_contours, motif_contours):
    #populates self.scores[alt] with contour values for the current residue
    #these contour values are: cablam, c_alpha_geom, alpha, beta, threeten
    category = self.contour_category()
    cablam_point = [self.measures. mu_in,self.measures.mu_out, self.measures.nu]
    ca_point = [self.measures.mu_in, self.measures.mu_out, self.measures.ca_virtual]
    motif_point = [self.measures.mu_in, self.measures.mu_out]
    if None in cablam_point: cablam=None
    else: cablam = cablam_contours[category].valueAt(cablam_point)
    if None in ca_point: c_alpha_geom=None
    else: c_alpha_geom = ca_contours[category].valueAt(ca_point)
    if None in motif_point: alpha, beta, threeten = 0, 0, 0
    else:
      alpha =    motif_contours['loose_alpha'].valueAt(motif_point)
      beta =     motif_contours['regular_beta'].valueAt(motif_point)
      threeten = motif_contours['loose_threeten'].valueAt(motif_point)
    self.scores = cablam_score(
      cablam=cablam,
      c_alpha_geom=c_alpha_geom,
      alpha=alpha,
      beta=beta,
      threeten=threeten)
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ set_cablam_feedback
  #-----------------------------------------------------------------------------
  def set_cablam_feedback(self):
    #populates self.feedback with True/False outlier status for easy access
    #these statuses are: cablam_outlier, cablam_disfavored, c_alpha_geom_outlier
    self.feedback = cablam_feedback()
    #outlier flags
    if self.scores.cablam is not None and self.scores.cablam < CABLAM_OUTLIER_CUTOFF:
      self.feedback.cablam_outlier = True
    else:
      self.feedback.cablam_outlier = False
    if self.scores.cablam is not None and self.scores.cablam < CABLAM_DISFAVORED_CUTOFF:
      self.feedback.cablam_disfavored = True
    else:
      self.feedback.cablam_disfavored = False
    if self.scores.c_alpha_geom is not None and self.scores.c_alpha_geom < CA_GEOM_CUTOFF:
      self.feedback.c_alpha_geom_outlier = True
    else:
      self.feedback.c_alpha_geom_outlier = False
    if self.feedback.cablam_outlier or self.feedback.cablam_disfavored or self.feedback.c_alpha_geom_outlier:
      self.outlier = True
    #secondary structure
    #This is semi-duplicated from assemble_secondary_structure
    if not self.prevres or not self.nextres:
      #alpha, beta, and threeten defaults are already None
      return
    if ((self.scores.alpha >= ALPHA_CUTOFF or self.scores.threeten >=THREETEN_CUTOFF)
      and (self.prevres.scores.alpha >= ALPHA_CUTOFF or self.prevres.scores.threeten >= THREETEN_CUTOFF)
      and (self.nextres.scores.alpha >= ALPHA_CUTOFF or self.nextres.scores.threeten >= THREETEN_CUTOFF)):
      if (self.scores.threeten > self.scores.alpha and
        (self.nextres.scores.threeten > self.nextres.scores.alpha or
          self.prevres.scores.threeten > self.prevres.scores.alpha)):
        self.feedback.threeten=True
      else:
        self.feedback.alpha=True
    if self.scores.beta >= BETA_CUTOFF and self.prevres.scores.beta >= BETA_CUTOFF and self.nextres.scores.beta >= BETA_CUTOFF:
      self.feedback.beta=True
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_kinemage
  #-----------------------------------------------------------------------------
  def as_kinemage(self, mode=None, out=sys.stdout):
    #prints kinemage markup for this residue
    #has separate output modes for cablam outliers and for ca geom outliers
    if mode == 'ca_geom':
      if self.feedback.c_alpha_geom_outlier is not None:
        stats = self.mp_id() + " ca_geom=%.2f alpha=%.2f beta=%.2f three-ten=%.2f" %(self.scores.c_alpha_geom*100, self.scores.alpha*100, self.scores.beta*100, self.scores.threeten*100)
        CA_1 = self.prevres.get_atom(' CA ').xyz
        CA_2 = self.get_atom(' CA ').xyz
        CA_3 = self.nextres.get_atom(' CA ').xyz
        out.write('\n{'+stats+'} P '+str(CA_2[0]-(CA_2[0]-CA_1[0])*0.9)+' '+str(CA_2[1]-(CA_2[1]-CA_1[1])*0.9)+' '+str(CA_2[2]-(CA_2[2]-CA_1[2])*0.9))
        out.write('\n{'+stats+'} '+str(CA_2[0])+' '+str(CA_2[1])+' '+str(CA_2[2]))
        out.write('\n{'+stats+'} '+str(CA_2[0]-(CA_2[0]-CA_3[0])*0.9)+' '+str(CA_2[1]-(CA_2[1]-CA_3[1])*0.9)+' '+str(CA_2[2]-(CA_2[2]-CA_3[2])*0.9))
    elif mode == 'cablam':
      if self.feedback.cablam_outlier is not None:
        stats = self.mp_id() + " cablam=%.2f alpha=%.2f beta=%.2f three-ten=%.2f" %(self.scores.cablam*100, self.scores.alpha*100, self.scores.beta*100, self.scores.threeten*100)
        CA_1, O_1 = self.prevres.get_atom(' CA ').xyz,self.prevres.get_atom(' O  ').xyz
        CA_2, O_2 = self.get_atom(' CA ').xyz,self.get_atom(' O  ').xyz
        CA_3      = self.nextres.get_atom(' CA ').xyz
        X_1 = perptersect(CA_1,CA_2,O_1)
        X_2 = perptersect(CA_2,CA_3,O_2)
        midpoint = [ (X_1[0]+X_2[0])/2.0 , (X_1[1]+X_2[1])/2.0 , (X_1[2]+X_2[2])/2.0 ]
        out.write('\n{'+stats+'} P '+ str(O_1[0]) +' '+ str(O_1[1]) +' '+ str(O_1[2]))
        out.write('\n{'+stats+'} '+ str(X_1[0]) +' '+ str(X_1[1]) +' '+ str(X_1[2]))
        out.write('\n{'+stats+'} '+ str(midpoint[0]) +' '+ str(midpoint[1]) +' '+ str(midpoint[2]))
        out.write('\n{'+stats+'} '+ str(X_2[0]) +' '+ str(X_2[1]) +' '+ str(X_2[2]))
        out.write('\n{'+stats+'} '+ str(O_2[0]) +' '+ str(O_2[1]) +' '+ str(O_2[2]))
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ is_same_as_other_result
  #-----------------------------------------------------------------------------
  def is_same_as_other_result(self,other_result):
    #Compare this result object to another to see if they are effectively the
    #  same.  Identical cablam geometry (mu_in, mu_out, and nu) is assumed to
    #  mean identical residues:
    if (self.measures.mu_in != other_result.measures.mu_in
      or self.measures.mu_out != other_result.measures.mu_out
      or self.measures.nu != other_result.measures.nu):
      return False
    return True
  #-----------------------------------------------------------------------------
  #}}}
#-------------------------------------------------------------------------------
#}}}

#{{{ cablamalyze class
#-------------------------------------------------------------------------------
class cablamalyze(validation):
  """
  Frontend for calculating cablam statistics for a model
  """
  __slots__ = validation.__slots__ + [
    "residue_count",
    "outlier_count",
    "sec_struc_records",
    "out",
    "pdb_hierarchy",
    "all_results"
    ]

  def get_result_class(self): return cablam_result

  #{{{ __init__
  #-----------------------------------------------------------------------------
  def __init__(self,
    pdb_hierarchy,
      outliers_only,
      out,
      quiet):
    validation.__init__(self)
    from mmtbx.validation import utils
    from scitbx.array_family import flex
    #self._outlier_i_seqs = flex.size_t()
    self.out = out
    cablam_contours = fetch_peptide_expectations()
    ca_contours = fetch_ca_expectations()
    motif_contours = fetch_motif_contours()
    self.pdb_hierarchy = pdb_hierarchy
    pdb_atoms = pdb_hierarchy.atoms()
    all_i_seqs = pdb_atoms.extract_i_seq()
    if all_i_seqs.all_eq(0):
      pdb_atoms.reset_i_seq()
    use_segids = utils.use_segids_in_place_of_chainids(
      hierarchy=pdb_hierarchy)

    self.all_results = {}
    ordered_keys = {}
    all_keys = []
    confs = []
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        if not chain.is_protein():
          continue
        for conf in chain.conformers():
          if conf.is_protein(): break #at least one conformer must be protein
        else: continue
        if use_segids:
          chain_id = utils.get_segid_as_chainid(chain=chain).rjust(2)
        else:
          chain_id = chain.id.rjust(2)
        #The above .rjust(2)'s are to force 2-char chain ids
        current_chain = cablam_chain()
        self.all_results[chain_id] = current_chain
        previous_result = None
        for conf in chain.conformers():
          if not conf.is_protein():
            continue
          current_conf = cablam_conf()
          current_conf.conf_name = conf.altloc
          current_chain.confs[conf.altloc] = current_conf
          current_chain.conf_names.append(conf.altloc)
          for residue in conf.residues():
            result = cablam_result(
              residue=residue,
              resseq=residue.resseq,
              icode=residue.icode,
              alt=conf.altloc,
              chain=chain_id,
              #formatting note: residue.id_str() = 'pdbres="THR B 182 "'
              #chain=residue.id_str()[11:13],
              #residue.id_str() turned out to break on some segid formatting
              prevres=None,
              nextres=None,
              has_ca=False,
              has_mc=False,
              outlier = False,
              measures=None,
              scores=None
              )
            result.check_atoms()
            result.link_residues(previous_result)
            #Occasionally a conformer may have more than one "residue" that has
            # the same sorting_id (sorting_id is just chain+resseq+icode)
            # see phenix_regression/pdb/lysozyme_nohoh_plus6H.pdb, where WAT 14
            # and ARG 14 have the same sorting_id
            #This check my not be the correct behavior to catch such a
            # formatting error.
            if result.sorting_id() not in current_conf.results:
              current_conf.results[result.sorting_id()] = result
            previous_result = result
    for chain in self.all_results:
      for conf in self.all_results[chain].confs:
        for result in self.all_results[chain].confs[conf].results:
          if self.all_results[chain].confs[conf].results[result].has_ca:
            self.all_results[chain].confs[conf].results[result].calculate_cablam_geometry()
            self.all_results[chain].confs[conf].results[result].calculate_contour_values(cablam_contours, ca_contours, motif_contours)
    for chain in self.all_results:
      for conf in self.all_results[chain].confs:
        for result in self.all_results[chain].confs[conf].results:
          if self.all_results[chain].confs[conf].results[result].has_ca:
            self.all_results[chain].confs[conf].results[result].set_cablam_feedback()
    self.assemble_secondary_structure()
    self.make_single_results_object(confs, all_keys)
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ assemble_secondary_structure
  #-----------------------------------------------------------------------------
  def assemble_secondary_structure(self):
    #assembles complete secondary structure elements (alpha helices, three-ten
    #  helices, and beta strands) from individual residue
    #have to assemble from scores to handle helix transitions
    for chain in self.all_results:
      for conf_id in self.all_results[chain].confs:
        conf = self.all_results[chain].confs[conf_id]
        records = []
        record_start = None
        helix_in_progress = False
        result_ids = conf.results.keys()
        result_ids.sort()
        for result_id in result_ids:
          result = conf.results[result_id]
          #is it evaluable?
          if not result.prevres:
            continue
          if not result.has_ca or not result.nextres:
            if helix_in_progress:
              records.append(secondary_structure_segment(start=record_start, end=result.prevres, segment_type=helix_in_progress, segment_length=record_length))
              helix_in_progress = False
            continue

          #helix building
          #This requires that the residues be the center of three in any combination of helix types
          #threeten segments of only 2 are lost in this method relative to the previous
          if ((result.scores.alpha >= ALPHA_CUTOFF or result.scores.threeten >=THREETEN_CUTOFF)
            and (result.prevres.scores.alpha >= ALPHA_CUTOFF or result.prevres.scores.threeten >= THREETEN_CUTOFF)
            and (result.nextres.scores.alpha >= ALPHA_CUTOFF or result.nextres.scores.threeten >= THREETEN_CUTOFF)):
            #now determine which helix type the current residue should be identified as
            #if at least two residues in a row have higher threeten scores than alpha scores, they can be considered threeten
            if (result.scores.threeten > result.scores.alpha and
              (result.nextres.scores.threeten > result.nextres.scores.alpha or
                result.prevres.scores.threeten > result.prevres.scores.alpha)):
              thisres = 'threeten'
            else:
              thisres = 'alpha'
            if helix_in_progress:
              if thisres == helix_in_progress: #is it same as previous residue
                record_length += 1
                continue
              else: #or has it changed helix types
                records.append(secondary_structure_segment(start=record_start, end=result.prevres, segment_type=helix_in_progress, segment_length=record_length))
                helix_in_progress = thisres
                record_start = result
                record_length = 1
            else:
              helix_in_progress = thisres
              record_start = result
              record_length = 1
          else: #(current residue is not helix)
            if helix_in_progress:
              #might fail on chain breaks?
              records.append(secondary_structure_segment(start=record_start, end=result.prevres, segment_type=helix_in_progress, segment_length=record_length))
              helix_in_progress = False
              record_start = None
              record_length = 0
            else:
              continue
        #helix building end

        #beta strands require another, separate pass
        strand_in_progress = False
        record_start = None
        for result_id in result_ids:
          result = conf.results[result_id]
          if not result.prevres:
            continue
          if not result.has_ca or not result.nextres:
            if strand_in_progress:
              records.append(secondary_structure_segment(start=record_start, end=result.prevres, segment_type='beta', segment_length=record_length))
              strand_in_progress = False
            continue
          if result.scores.beta >= BETA_CUTOFF and result.prevres.scores.beta >= BETA_CUTOFF and result.nextres.scores.beta >= BETA_CUTOFF:
            if strand_in_progress:
              record_length += 1
              continue
            else:
              strand_in_progress = True
              record_start = result
              record_length = 1
          else:
            if strand_in_progress:
              records.append(secondary_structure_segment(start=record_start, end=result.prevres, segment_type='beta', segment_length=record_length))
              strand_in_progress = False
              record_start = None
              record_length = 0
            else:
              continue
        #beta strand building end
        #NOTE: Each strand is currently treated as its own sheet
        #Developing or implementing strand-to-sheet assembly that does not rely on
        #  H-bonds is a major future goal
        conf.sec_struc_records = records
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ make_single_results_object
  #-----------------------------------------------------------------------------
  def make_single_results_object(self, confs, all_keys):
    #should work without any arguments
    #populates self.results
    self.results = []
    chains = self.all_results.keys()
    chains.sort()
    for chain_id in chains:
      chain = self.all_results[chain_id]
      #take the first conformer as the basis for comparison
      conf  = chain.confs[chain.conf_names[0]]
      if len(chain.conf_names) == 0:
        for result_id in conf.results:
          result = conf.results[result_id]
          result.alt = ''
          #set self.results id
        continue #go to next chain
      #else, combine results into single list
      result_ids = conf.results.keys()
      result_ids.sort()
      #for result_id in conf.results:
      for result_id in result_ids:
        result = conf.results[result_id]
        if not result.has_ca: continue
        #results without CAs have measures=None and break the
        #  is_same_as_other_result check. Also, they aren't evaluable residues.
        self.results.append(result)
        found_meaningful_alt = False
        for other_conf in chain.conf_names[1:]:
          if result.sorting_id() in chain.confs[other_conf].results:
            other_result = chain.confs[other_conf].results[result.sorting_id()]
            if not other_result.has_ca: continue
            if not result.is_same_as_other_result(other_result):
              self.results.append(other_result)
              found_meaningful_alt = True
        if not found_meaningful_alt:
          result.alt = ''
          #set self.results id
          pass
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_text
  #-----------------------------------------------------------------------------
  def as_text(self, outliers_only=False):
    #prints colon-separated text for CaBLAM validation
    #one line per residue (or alternate)
    #Output is formatted to be human-readable, and is also used by MolProbity
    #This is the default output for running this script from commandline
    self.out.write('residue : outlier_type : contour_level : ca_contour_level : sec struc recommendation : alpha score : beta score : three-ten score')
    for result in self.results:
      if not result.has_ca:
        continue
      #if not result.feedback:
      #  continue
      if outliers_only:
        if not (result.feedback.cablam_disfavored or result.feedback.c_alpha_geom_outlier):
          continue

      if result.feedback.c_alpha_geom_outlier:
        outlier_type = ' CA Geom Outlier    '
      elif result.feedback.cablam_outlier:
        outlier_type = ' CaBLAM Outlier     '
      elif result.feedback.cablam_disfavored:
        outlier_type = ' CaBLAM Disfavored  '
      else:
        outlier_type = '                    '

      if result.scores.cablam is not None:
        cablam_level = '%.5f' %result.scores.cablam
      else:
        cablam_level = '       ' #default printing for CA-only models
      if result.scores.c_alpha_geom is not None:
        ca_geom_level = '%.5f' %result.scores.c_alpha_geom
      else:
        continue #if this is missing, there's nothing

      if result.feedback.c_alpha_geom_outlier:
        suggestion = '                 '
      elif result.feedback.threeten:
        suggestion = ' try three-ten   '
      elif result.feedback.alpha:
        suggestion = ' try alpha helix '
      elif result.feedback.beta:
        suggestion = ' try beta sheet  '
      else:
        suggestion = '                 '

      outlist = [result.mp_id() ,outlier_type, cablam_level, ca_geom_level, suggestion, '%.5f' %result.scores.alpha, '%.5f' %result.scores.beta, '%.5f' %result.scores.threeten]
      #outlist = [result.residue.id_str() ,outlier_type, cablam_level, ca_geom_level, suggestion, '%.5f' %result.scores.alpha, '%.5f' %result.scores.beta, '%.5f' %result.scores.threeten]
      self.out.write('\n'+':'.join(outlist))
    self.out.write('\n')
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_oneline
  #-----------------------------------------------------------------------------
  def as_oneline(self):
    #prints a one-line summary of cablam statistics for a structure
    #for oneline purposes, alternates are collapsed: each residue contributes up
    #  to 1 to each outlier count, regarless of how many outlier alternates it
    #  may contain
    residue_count = 0
    cablam_outliers = 0
    cablam_disfavored = 0
    ca_geom_outliers = 0
    prev_result_id = None
    is_cablam_outlier = 0
    is_cablam_disfavored = 0
    is_ca_geom_outlier = 0
    for result in self.results:
      if not result.has_ca:
        continue
      if result.scores.cablam is None:
        is_residue = 0
        continue
      is_residue = 1
      if result.sorting_id() != prev_result_id:
        #new residue; update counts
        residue_count    += is_residue
        cablam_outliers  += is_cablam_outlier
        cablam_disfavored+= is_cablam_disfavored
        ca_geom_outliers += is_ca_geom_outlier
        is_cablam_outlier    = 0
        is_cablam_disfavored = 0
        is_ca_geom_outlier   = 0
      if result.scores.cablam < CABLAM_OUTLIER_CUTOFF:
        is_cablam_outlier    = 1
      if result.scores.cablam < CABLAM_DISFAVORED_CUTOFF:
        is_cablam_disfavored = 1
      if result.scores.c_alpha_geom is not None and result.scores.c_alpha_geom < CA_GEOM_CUTOFF:
        is_ca_geom_outlier = 1
      prev_result_id = result.sorting_id()
    residue_count    += is_residue
    cablam_outliers  += is_cablam_outlier
    cablam_disfavored+= is_cablam_disfavored
    ca_geom_outliers += is_ca_geom_outlier
    if residue_count == 0:
      self.out.write('pdbid:0:0:0:0\n')
    else:
      self.out.write('pdbid:%i:%.1f:%.1f:%.2f\n' %(residue_count, cablam_outliers/residue_count*100, cablam_disfavored/residue_count*100, ca_geom_outliers/residue_count*100) )
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_kinemage
  #-----------------------------------------------------------------------------
  def as_kinemage(self):
    #output cablam validation as standalone kinemage markup for viewing in KiNG
    self.out.write('\n@subgroup {cablam disfavored} dominant\n')
    self.out.write('@vectorlist {cablam disfavored} color= purple width= 4 master={cablam disfavored} off') #default off
    for result in self.results:
      if not result.has_ca:
        continue
      if result.feedback.cablam_disfavored:
        result.as_kinemage(mode="cablam", out=self.out)
    self.out.write('\n@subgroup {cablam outlier} dominant\n')
    self.out.write('@vectorlist {cablam outlier} color= magenta width= 4 master={cablam outlier}') #default on
    for result in self.results:
      if not result.has_ca:
        continue
      if result.feedback.cablam_outlier:
        result.as_kinemage(mode="cablam", out=self.out)
    self.out.write('\n@subgroup {ca geom outlier} dominant\n')
    self.out.write('@vectorlist {ca geom outlier} color= red width= 4 master={ca geom outlier}') #default on
    for result in self.results:
      if not result.has_ca:
        continue
      if result.feedback.c_alpha_geom_outlier:
        result.as_kinemage(mode="ca_geom", out=self.out)
    self.out.write('\n')
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_full_kinemage
  #-----------------------------------------------------------------------------
  def as_full_kinemage(self,pdbid=''):
    #output cablam validation as kinemage markup on pdb model. Future version
    #  of this will also print ribbons, but standalone ribbon code must be
    #  developed first.
    #Pdb-to-kinemage printing has been hijacked from mmtbx.kinemage.validation
    #  That code was not meant to be used outside its original context, so this
    #  may be fragile.
    from mmtbx.kinemage.validation import get_kin_lots, build_name_hash
    from mmtbx import monomer_library
    from mmtbx.monomer_library import pdb_interpretation
    i_seq_name_hash = build_name_hash(pdb_hierarchy=self.pdb_hierarchy)
    sites_cart=self.pdb_hierarchy.atoms().extract_xyz()
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
    pdb_io=self.pdb_hierarchy.as_pdb_input()
    processed_pdb_file = pdb_interpretation.process(
        mon_lib_srv=mon_lib_srv,
        ener_lib=ener_lib,
        pdb_inp=pdb_io,
        #params=work_params.kinemage.pdb_interpretation,
        substitute_non_crystallographic_unit_cell_if_necessary=True)
    geometry = processed_pdb_file.geometry_restraints_manager()
    flags = geometry_restraints.flags.flags(default=True)
    #angle_proxies = geometry.angle_proxies
    pair_proxies = geometry.pair_proxies(flags=flags, sites_cart=sites_cart)
    bond_proxies = pair_proxies.bond_proxies
    quick_bond_hash = {}
    for bp in bond_proxies.simple:
      if (i_seq_name_hash[bp.i_seqs[0]][9:14] == i_seq_name_hash[bp.i_seqs[1]][9:14]):
        if quick_bond_hash.get(bp.i_seqs[0]) is None:
          quick_bond_hash[bp.i_seqs[0]] = []
        quick_bond_hash[bp.i_seqs[0]].append(bp.i_seqs[1])

    for model in self.pdb_hierarchy.models():
      for chain in model.chains():
        self.out.write(get_kin_lots(chain, bond_hash=quick_bond_hash, i_seq_name_hash=i_seq_name_hash, pdbID=pdbid, index=0, show_hydrogen=True))
    self.as_kinemage()
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_secondary_structure
  #-----------------------------------------------------------------------------
  def as_secondary_structure(self, conf_request=None):
    #returns CaBLAM secondary structure identification as phenix's preferred
    #  iotbx.secondary_structure objects
    #each individual beta strand is currently represented as a whole sheet
    #Proper sheet reporting will depend on developing or implementing a
    #  strand-to-sheet assembly that does not require H-bonds
    from iotbx.pdb import secondary_structure
    helix_i = 0
    sheet_i = 0
    helix_records = []
    strand_records = []
    #return_records = []

    chain_list = self.all_results.keys()
    chain_list.sort()
    for chain_id in chain_list:
      chain = self.all_results[chain_id]
      if conf_request in chain.conf_names:
        conf = conf_request
      else:
        conf = chain.conf_names[0]
      for record in chain.confs[conf].sec_struc_records:
        if record.segment_type == 'alpha' or record.segment_type == 'threeten':
          helix_i += 1
          if record.segment_type == 'alpha':
            helix_class = 1
          elif record.segment_type == 'threeten':
            helix_class = 5
          return_record = secondary_structure.pdb_helix(
            serial = helix_i,
            helix_id = helix_i,
            start_resname  = record.start.residue.resname,
            start_chain_id = record.start.chain,
            #start_chain_id = " A",
            start_resseq   = record.start.resseq,
            start_icode    = record.start.icode,
            end_resname    = record.end.residue.resname,
            end_chain_id   = record.end.chain,
            end_resseq     = record.end.resseq,
            end_icode      = record.end.icode,
            helix_class    = helix_class,
            comment = "",
            length = record.segment_length)
          helix_records.append(return_record)
          #return_records.append(return_record)
        if record.segment_type == 'beta':
          sheet_i += 1
          strand_record = secondary_structure.pdb_strand(
            sheet_id       = sheet_i,
            strand_id      = 1,
            start_resname  = record.start.residue.resname,
            start_chain_id = record.start.chain,
            start_resseq   = record.start.resseq,
            start_icode    = record.start.icode,
            end_resname    = record.end.residue.resname,
            end_chain_id   = record.end.chain,
            end_resseq     = record.end.resseq,
            end_icode      = record.end.icode,
            sense          = 1
            )
          return_record=secondary_structure.pdb_sheet(
            sheet_id = sheet_i,
            n_strands = 1,
            strands = [strand_record],
            registrations = [None],
            hbond_list = []
            )
          strand_records.append(return_record)
          #return_records.append(return_record)
    return secondary_structure.annotation(helices=helix_records,sheets=strand_records)
    #return return_records
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_records
  #-----------------------------------------------------------------------------
  def as_records(self, conf_request=None):
    #outputs pdb-style HELIX and SHEET secondary structure records
    #uses the iotbx.secondary_structure object
    #By default, this returns the first conformation (alt) in each chain
    #Other conformations can be accessed with conf_request
    records = self.as_secondary_structure(conf_request=conf_request)
    self.out.write(records.as_pdb_str()+'\n')
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_records_and_pdb
  #-----------------------------------------------------------------------------
  def as_records_and_pdb(self, conf_request=None):
    #outputs pdb-style HELIX and SHEET secondary structure records, followed by
    #  the pdb file
    self.as_records(conf_request=conf_request)
    self.out.write(self.pdb_hierarchy.as_pdb_string())
  #-----------------------------------------------------------------------------
  #}}}
#-------------------------------------------------------------------------------
#}}}

#{{{ run
#-------------------------------------------------------------------------------
def run(args):
  #{{{ phil parsing
  #-----------------------------------------------------------------------------
  interpreter = libtbx.phil.command_line.argument_interpreter(master_phil=master_phil)
  sources = []
  for arg in args:
    if os.path.isfile(arg): #Handles loose filenames
      input_file = file_reader.any_file(arg)
      if (input_file.file_type == "pdb"):
        sources.append(interpreter.process(arg="pdb_infile=\"%s\"" % arg))
      elif (input_file.file_type == "phil"):
        sources.append(input_file.file_object)
    else: #Handles arguments with xxx=yyy formatting
      arg_phil = interpreter.process(arg=arg)
      sources.append(arg_phil)
  work_phil = master_phil.fetch(sources=sources)
  work_params = work_phil.extract()
  params = work_params.cablam
  #-----------------------------------------------------------------------------
  #}}}

  if params.help:
    usage()
    interpretation()
    sys.exit()

  if not params.pdb_infile:
    sys.stderr.write(
      '\nMissing input data, please provide .pdb file\n')
    usage()
    sys.exit()

  if not os.path.isfile(params.pdb_infile):
    sys.stderr.write(params.pdb_infile + " is not a file or could not be found")
    sys.exit()
  else:
    pdb_infile = params.pdb_infile

  pdb_in = file_reader.any_file(pdb_infile)
  if pdb_in.file_type != "pdb":
    sys.stderr.write(pdb_infile +" not id'd as readable file\n")
    sys.exit()
  pdbid = os.path.splitext(os.path.basename(pdb_infile))[0]
  pdb_io = pdb.input(pdb_infile)
  input_hierarchy = pdb_io.construct_hierarchy()
  cablam = cablamalyze(
    pdb_hierarchy=input_hierarchy,
    outliers_only=False,
    out=sys.stdout,
    quiet=False)

  #output = *text kin full_kin records records_and_pdb oneline
  if params.output=='oneline':
    cablam.as_oneline()
  elif params.output=='kin':
    cablam.as_kinemage()
  elif params.output=='full_kin':
    cablam.as_full_kinemage(pdbid=pdbid)
  elif params.output=='records':
    cablam.as_records()
  elif params.output=='records_and_pdb':
    cablam.as_records_and_pdb()
  else: #default text output
    cablam.as_text(outliers_only=params.outliers_only)
#-------------------------------------------------------------------------------
#}}}

if __name__ == "__main__":
  run(sys.argv[1:])
