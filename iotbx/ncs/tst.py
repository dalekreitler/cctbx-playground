import iotbx.ncs
import iotbx.pdb
import iotbx.phil
from libtbx.test_utils import approx_equal

def exercise_00(prefix="iotbx_ncs_exercise_00"):
  pdb_file_name = "%s.pdb"%prefix
  pdb_str_1 = """
CRYST1  399.000  399.000  399.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A  34     125.208 211.886 175.417  1.00  0.00           N
ATOM      2  CA  GLY A  34     125.035 211.123 174.168  1.00  0.00           C
ATOM      3  C   GLY A  34     126.386 210.806 173.507  1.00  0.00           C
ATOM      4  O   GLY A  34     127.304 211.628 173.503  1.00  0.00           O
TER
ATOM      5  N   GLY B  34     251.532 143.432 175.422  1.00  0.00           N
ATOM      6  CA  GLY B  34     252.120 143.948 174.173  1.00  0.00           C
ATOM      7  C   GLY B  34     251.212 144.998 173.512  1.00  0.00           C
ATOM      8  O   GLY B  34     249.986 144.872 173.510  1.00  0.00           O
TER
ATOM      9  N   GLY C  34     189.583 273.076 175.423  1.00  0.00           N
ATOM     10  CA  GLY C  34     188.804 273.006 174.173  1.00  0.00           C
ATOM     11  C   GLY C  34     188.920 271.622 173.510  1.00  0.00           C
ATOM     12  O   GLY C  34     189.986 271.004 173.508  1.00  0.00           O
TER
  """
  pdb_str_2 = """
CRYST1  399.000  399.000  399.000  90.00  90.00  90.00 P 1
ATOM      1  O   HOH S   1     109.583 203.076 175.423  1.00  0.00           O
TER
ATOM      1  N   GLY A  34     125.208 211.886 175.417  1.00  0.00           N
ATOM      2  CA  GLY A  34     125.035 211.123 174.168  1.00  0.00           C
ATOM      3  C   GLY A  34     126.386 210.806 173.507  1.00  0.00           C
ATOM      4  O   GLY A  34     127.304 211.628 173.503  1.00  0.00           O
TER
ATOM      5  N   GLY B  34     251.532 143.432 175.422  1.00  0.00           N
ATOM      6  CA  GLY B  34     252.120 143.948 174.173  1.00  0.00           C
ATOM      7  C   GLY B  34     251.212 144.998 173.512  1.00  0.00           C
ATOM      8  O   GLY B  34     249.986 144.872 173.510  1.00  0.00           O
TER
ATOM      9  N   GLY C  34     189.583 273.076 175.423  1.00  0.00           N
ATOM     10  CA  GLY C  34     188.804 273.006 174.173  1.00  0.00           C
ATOM     11  C   GLY C  34     188.920 271.622 173.510  1.00  0.00           C
ATOM     12  O   GLY C  34     189.986 271.004 173.508  1.00  0.00           O
TER
ATOM      9  O   TYR D   4     189.583 273.076 175.423  1.00  0.00           O
TER
  """
  ncs_params_str = """
ncs_group_selection {
  ncs_group {
    master_ncs_selection = chain A
    selection_copy = chain B
    selection_copy = chain C
  }
}
  """
  def check_result(ncs_inp):
    ncs_groups = ncs_inp.get_ncs_restraints_group_list()
    assert len(ncs_groups) == 1
    ncs_group = ncs_groups[0]
    assert approx_equal(ncs_group.master_ncs_iselection, [0,1,2,3])
    assert len(ncs_group.copies) == 2
    assert approx_equal(ncs_group.copies[0].ncs_copy_iselection, [4,5,6,7])
    assert approx_equal(ncs_group.copies[1].ncs_copy_iselection, [8,9,10,11])
  for pdb_str in [pdb_str_1, pdb_str_2]:
    of = open(pdb_file_name, "w")
    print >> of, pdb_str
    of.close()
    pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
    # using pdb_inp
    ncs_inp = iotbx.ncs.input(pdb_inp = pdb_inp)
    check_result(ncs_inp)
    # using file_name
    ncs_inp = iotbx.ncs.input(file_name = pdb_file_name)
    check_result(ncs_inp)
    # using pdb string
    ncs_inp = iotbx.ncs.input(pdb_string = pdb_str)
    check_result(ncs_inp)
    # using combination of pdb_inp and Phil parameter string
    ncs_inp = iotbx.ncs.input(pdb_inp = pdb_inp,
      ncs_selection_params = ncs_params_str)
    check_result(ncs_inp)
    # using combination of pdb file name and Phil parameter string
    ncs_inp = iotbx.ncs.input(file_name = pdb_file_name,
      ncs_selection_params = ncs_params_str)
    check_result(ncs_inp)
    # using combination of pdb string and Phil parameter string
    ncs_inp = iotbx.ncs.input(pdb_string = pdb_str,
      ncs_selection_params = ncs_params_str)
    check_result(ncs_inp)

if (__name__ == "__main__"):
  exercise_00()
