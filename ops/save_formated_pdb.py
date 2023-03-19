from pymol import cmd
import sys
import os

#pymol -cq xx.py [pdb_path]
mobile_name = sys.argv[3]
save_path = sys.argv[4]
mobile_name = os.path.abspath(mobile_name)
mobile_state=os.path.split(mobile_name)[1][:-4]
cmd.load(mobile_name,mobile_state)
#save_path = mobile_name[:-4]+"_formated.pdb"
cmd.save(save_path,mobile_state)
