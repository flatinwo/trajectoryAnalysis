#!/usr/bin/env python

import numpy as np
import os

def molecularize(mol,size=4):
	""" a simple function to compute center of mass
	of a given molecule

	params: mol should be a list"""

	assert size == len(mol)
	# split molecule
	s_molecule = np.array([ [float(i) for i in atom.split()[1:]] for atom in mol ])
	return (mol[0][0],s_molecule.mean(axis=0))

def splitxyzfile(fn,molsize=4):
	with open(fn) as f:
		content = f.read().splitlines()

	header, body = content[0:2],content[2:]
	molecules = list(map(list,zip(*[iter(body)]*molsize)))
	return header, molecules

def get_min_idx(molecularized):
	mmin = np.inf
	min_idx = np.inf
	for i in range(len(molecularized)):
		current_min = np.linalg.norm(molecularized[i][1])
		if current_min < mmin : 
			min_idx = i
			mmin = current_min 
	return min_idx

def processbop(fn,ntype=2,molsize=4,folders=['self','cross','com']):
	header, molecules = splitxyzfile(fn)
	type_i_molecules = [ [ molecularize(molecule) for molecule in molecules[i::ntype] ] for i in range(0,ntype) ]
	filenames = [ 'molecules_'+str(i)+'.xyz' for i in range(0,len(type_i_molecules))]

	#make folders for calculations
	for f in folders: os.makedirs(f,exist_ok=True)

	#write molecules
	i=0
	for fn in filenames:
		header[0]=str(molsize*len(molecules[i::ntype]))
		with open('/'.join(['.',folders[0],fn]),'w') as f:
			f.writelines('{}\n'.format('\n'.join(header)))
			f.writelines('{}\n'.format('\n'.join(line)) for line in molecules[i::ntype])
		i+=1
			
	min_indices = []
	for mtype,fn in zip(type_i_molecules,filenames):
		min_idx = get_min_idx(mtype)
		lines =  [[str(j) for j in sum([ list(i) for i in atom ],[]) ] for atom in mtype]
		header[0] = str(len(lines))

		min_indices.append(min_idx)

		#write center of mass
		with open('/'.join(['.',folders[2],'com_'+fn]),'w') as f:
			f.writelines('{}\n'.format('\n'.join(header)))
			f.writelines('{}\n'.format('\t'.join(line)) for line in lines)

	
	# write cross molecules
	# brute force
	i=1
	print(min_indices)
	import itertools as IT
	for fn,idx in zip(filenames,min_indices):
		header[0]=str(molsize*len(molecules[(i+1)%2::ntype]))
		counter = IT.count(0)
		with open('/'.join(['.',folders[1],str(i)+"_"+fn]),'w') as f:
			f.writelines('{}\n'.format('\n'.join(header)))
			f.writelines('{}\n'.format('\n'.join(line)) for line in molecules[(i+1)%2::ntype] if idx != next(counter))
			f.writelines('{}\n'.format('\n'.join(molecules[i::ntype][min_indices[i]])))
		i-=1
	
		
def calculatebops(folders,filenames,ls=[2+2*i for i in range(6)]):
	import subprocess	
	curr_dir=os.getcwd()
	for f,fn in zip(folders,filenames):
		os.chdir('/'.join([curr_dir,f]))
		for l in ls: subprocess.run(["chiralityBOP", fn, str(l)])
  
		
if __name__ == "__main__":
	folders=['BOPs/self','BOPs/cross','BOPs/com']
	filenames = ['molecules_0.xyz','1_molecules_0.xyz']
	try:
		os.chdir(os.getcwd())
		processbop(fn="./best_lattices.xyz",folders=folders)
		calculatebops(folders=folders[0:2],filenames=filenames)
	except OSError:
		raise
		print("File not found")
	
	

