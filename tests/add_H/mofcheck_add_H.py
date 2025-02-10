from mofchecker import MOFChecker
import os

cif_path = './add_H_cifs'
cifs = os.listdir(cif_path)
cifs.sort()
for structure in cifs:
    mofchecker = MOFChecker.from_cif(os.path.join(cif_path,structure), h_number=4)
    mofchecker.adding_hydrogen
