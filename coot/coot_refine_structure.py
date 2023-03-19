import mrcfile
import numpy as np
import os
import shutil
def coot_refine_structure(input_cif_path,input_map_path,output_cif_path,coot_software="coot"):
    """
    input and output must be pdb format
    """
    output_cif_path = os.path.abspath(output_cif_path)
    work_dir = os.path.split(output_cif_path)[0]
    script_path= os.path.join(os.getcwd(),"coot")
    script_path = os.path.join(script_path,"all_atom_refine.py")
    new_script_path = os.path.join(work_dir,"all_atom_refine.py")
    shutil.copy(script_path,new_script_path)
    #determine contour
    with mrcfile.open(input_map_path,permissive=True) as mrc:
        data=mrc.data
    min_value= np.min(data[data>0])
    script_path = os.path.join(work_dir,"coot_refine.py")
    with open(script_path,'w') as file:
        #file.write('from all_atom_refine import real_space_atom_refine_output\n')
        with open(new_script_path,'r') as rfile:
            for line in rfile:
                file.write(line)
        file.write("\n")
        file.write('real_space_atom_residuerefine_output("%s","%s","%s",%f)'%(input_map_path,input_cif_path,output_cif_path,min_value))

    root_dir=os.getcwd()
    os.chdir(work_dir)
    os.system("%s --no-graphics --script coot_refine.py"%coot_software)
    os.chdir(root_dir)
