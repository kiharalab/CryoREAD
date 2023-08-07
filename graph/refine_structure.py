import os
from ops.os_operation import mkdir
import shutil
from data_processing.format_pdb import format_pdb,remove_op3_pdb
def refine_structure(input_pdb,input_map,output_dir,params):
    mkdir(output_dir)
    assert input_pdb[-4:]==".pdb"
    format_pdb_path =  os.path.join(output_dir,"input_format.pdb")
    format_pdb(input_pdb,format_pdb_path)
    if "resolution" not in params:
        params['resolution'] =5
        print("no resolutino input detected, use 5A as default resolution for refinement!")
    os.system('cd %s; phenix.real_space_refine %s %s resolution=%.4f '
              'output.suffix="_phenix_refine" '
              'skip_map_model_overlap_check=True'%(output_dir,
                                                   format_pdb_path,input_map,params['resolution']))
    gen_pdb_path = format_pdb_path[:-4]+"_phenix_refine_000.pdb"
    count_check=0
    while not os.path.exists(gen_pdb_path) and count_check<5:
        gen_pdb_path = format_pdb_path[:-4]+"_phenix_refine_00%d.pdb"%(count_check+1)
        count_check+=1
    if not os.path.exists(gen_pdb_path):
        print("1st round phenix refinement failed!")
        return
    refine1_pdb_path = os.path.join(output_dir,"Refine_cycle1.pdb")
    shutil.copy(gen_pdb_path,refine1_pdb_path)

    refine2_pdb_path = os.path.join(output_dir,"Refine_cycle2.pdb")
    from coot.coot_refine_structure import coot_refine_structure
    coot_software="coot"
    coot_refine_structure(refine1_pdb_path,input_map,refine2_pdb_path,coot_software)
    if not os.path.exists(refine2_pdb_path):
        print("2nd round coot refinement failed!")
        return

    refine3_pdb_path = os.path.join(output_dir,"Refine_cycle3.pdb")
    os.system('cd %s; phenix.real_space_refine %s %s resolution=%.4f '
                          'output.suffix="_phenix_refine" '
              'skip_map_model_overlap_check=True'%(output_dir,refine2_pdb_path,input_map,params['resolution']))
    phenix_final_pdb = refine2_pdb_path[:-4]+"_phenix_refine_000.pdb"
    count_check=0
    while not os.path.exists(phenix_final_pdb) and count_check<5:
        phenix_final_pdb = refine2_pdb_path[:-4]+"_phenix_refine_00%d.pdb"%(count_check+1)
        count_check+=1
    if not os.path.exists(phenix_final_pdb):
        print("3rd round phenix refinement failed!")
        return
    shutil.move(phenix_final_pdb,refine3_pdb_path)




