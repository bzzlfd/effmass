include ../makefile.include
all: upf2upfSO.x replace_ppbeta.x replace_ppchi.x

upf2upfSO.x: upf2upfSO.f90
	ifort -I$(IOTK_DIR) constant.f90 file_unit.f90 error_message.f90 \
		transform_to_upper.f90 read_key_words.f90 scan_key_words.f90 \
		gen_element_name_number.f90 upfpsp_module.f90 control_config.f90 atom_base.f90 \
		pse_pot_config.f90 \
		read_upf_v2.f90 cal_ityp.f90 check_atom_mov.f90 capital.f90 \
		read_vwr_head.f90 read_uspp_head.f90 read_atom_file.f90 \
		write_upf_v2.f90 upf2upfSO.f90 $(IOTK_DIR)libiotk.a \
		 -o $@

replace_ppbeta.x: main_replace_ppbeta.f90
	ifort -I$(IOTK_DIR) constant.f90 file_unit.f90 error_message.f90 \
		transform_to_upper.f90 read_key_words.f90 scan_key_words.f90 \
		gen_element_name_number.f90 upfpsp_module.f90 control_config.f90 atom_base.f90 \
		pse_pot_config.f90 \
		read_upf_v2.f90 cal_ityp.f90 check_atom_mov.f90 capital.f90 \
		read_vwr_head.f90 read_uspp_head.f90 read_atom_file.f90 \
		write_upf_v2.f90 main_replace_ppbeta.f90 $(IOTK_DIR)libiotk.a \
	        -o $@

replace_ppchi.x: main_replace_ppchi.f90
	ifort -I$(IOTK_DIR) constant.f90 file_unit.f90 error_message.f90 \
		transform_to_upper.f90 read_key_words.f90 scan_key_words.f90 \
		gen_element_name_number.f90 upfpsp_module.f90 control_config.f90 atom_base.f90 \
		pse_pot_config.f90  \
		read_upf_v2.f90 cal_ityp.f90 check_atom_mov.f90 capital.f90 \
		read_vwr_head.f90 read_uspp_head.f90 read_atom_file.f90 \
		write_upf_v2.f90 main_replace_ppchi.f90 $(IOTK_DIR)libiotk.a \
		-o $@



clean:
	rm -f *.x *.mod
