from mofchecker import MOFChecker
import os
import csv
import pandas as pd

def check(path):
	results = []
	descriptors = ['name','formula', 'spacegroup_symbol','spacegroup_number','symmetry_hash']

	output_file = 'check.csv'
	existing_files = []
	if os.path.exists(output_file):
	  df = pd.read_csv(output_file, header=None)
	  existing_files = set(df[0].astype(str).tolist())
	if not os.path.exists(output_file):
	  df = pd.DataFrame(columns=descriptors)
	  df.to_csv(output_file, index=False)
	for structure in path:
		str_name = structure.split('.cif')[0]
		if str_name not in existing_files:
		       mofchecker = MOFChecker.from_cif(os.path.join(cif_path,structure), primitive=False)
		       results = mofchecker.get_mof_descriptors(descriptors)
		       df = pd.DataFrame([results])
		       df.to_csv(output_file, mode='a', header=False, index=False)	

if __name__ == "__main__":
	cif_path = './cifs'
	cifs = os.listdir(cif_path)
	cifs.sort()
	check(cifs)   

