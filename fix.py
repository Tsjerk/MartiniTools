#!/usr/bin/env python

import sys
import modeller 
import modeller.scripts

env = modeller.environ()
env.io.atom_files_directory = ['../atom_files']
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

modeller.scripts.complete_pdb(env,sys.argv[1]).write(sys.stdout)


