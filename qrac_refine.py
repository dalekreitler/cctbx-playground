from iotbx import pdb
from iotbx.reflection_file_reader import any_reflection_file
from cctbx import xray
from cctbx.xray import structure_factors
from mmtbx.masks import masks
import mmtbx.f_model as fm
import argparse


command_line_parser = argparse.ArgumentParser()
command_line_parser.add_argument('pdb_input', action='store')
command_line_parser.add_argument('reflection_input', action='store')
command_line_args = command_line_parser.parse_args()


class refinement:


    def __init__(self,
		 refinement_type = 'None'):

        self.pdb_structure = pdb.hierarchy.input(command_line_args.pdb_input)
	self.get_crystal_symmetry()
	self.xray_structure = self.pdb_structure.input.xray_structure_simple()
	self.miller_arrays = any_reflection_file(command_line_args.reflection_input).as_miller_arrays()
	self.miller_arrays_work_test_partition()

        return
	
    def get_crystal_symmetry(self):
		
        if (self.pdb_structure):
	    self.refinement_symmetry = self.pdb_structure.input.crystal_symmetry()
	else:
	    self.generate_pdb_structure()
	    self.refinement_symmetry = self.pdb_structure.input.crystal_symmetry()
	return

    def miller_arrays_work_test_partition(self):

	'''generate miller array as amplitudes, 
	   separate work and test sets,
	   ensure working and test sets match dimensions,
	   expand to include Friedel mates'''

	f_obs = self.miller_arrays[0].map_to_asu()
	flags = self.miller_arrays[4]
	
	

	flags = flags.customized_copy(data=flags.data()=='f').map_to_asu()

	

	f_obs = f_obs.merge_equivalents().array()
	f_obs = f_obs.generate_bijvoet_mates()
	flags = flags.merge_equivalents().array()
	flags_plus_minus = flags.generate_bijvoet_mates()

	'''common_sets() ensures miller arrays are sorted,
	   are the same dimension, and have matching indices'''

	self.f_obs, self.flags_plus_minus = f_obs.common_sets(other=flags_plus_minus)
	self.f_obs_work = f_obs.select(~self.flags_plus_minus.data())
	self.f_obs_test = f_obs.select(self.flags_plus_minus.data())

	return

    def generate_f_calc(self,
                        miller_array = None,
                        xray_structure = None):

	xrs = xray_structure
	if (xrs is None): xrs = self.xray_structure
	if (miller_array is None): miller_array = self.f_obs
	manager = miller_array.structure_factors_from_scatterers(
            xray_structure = xrs,
            algorithm = 'fft',
            cos_sin_table = True) 
        
	self.f_calc = manager.f_calc()
	
	return

    def generate_f_model_manager(self):
		

        f_model_manager = fm.manager(f_calc = self.f_calc,
                                     f_obs = self.f_obs,
				     r_free_flags = self.flags_plus_minus,
				     xray_structure=self.xray_structure,
				     target_name='ls')

	f_model = f_model_manager.f_model()
	
	f_model_manager.update_solvent_and_scale()
	f_model_manager.show()
	
	k = f_model_manager.k_masks()
	print(dir(k))

	return


    def output_pdb_file(self):

        self.get_crystal_symmetry()
	output_pdb_file = open('test.pdb', 'w')
	output_pdb_file.write(self.pdb_structure.hierarchy.as_pdb_string(crystal_symmetry = self.refinement_symmetry))
	output_pdb_file.close()

	return

refine = refinement()
refine.generate_f_calc()
#refine.calculate_mask()
refine.generate_f_model_manager()

