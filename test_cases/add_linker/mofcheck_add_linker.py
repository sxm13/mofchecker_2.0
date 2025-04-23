from mofchecker import MOFChecker
import os

cif_path = './add_linker_cifs'
cifs = os.listdir(cif_path)
cifs.sort()
for structure in cifs:
    mofchecker = MOFChecker.from_cif(os.path.join(cif_path,structure),os.path.join('./linkers/NO3.cif'),linker_number=4)
    mofchecker.adding_linker
