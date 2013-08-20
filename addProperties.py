"""
    Module:  addProperties

    Reads a table and adds a column of the specified type.
    --------------------

    Felix Kruger
    momo.sander@googlemail.com
"""
import queryDevice
import os

def get_idx(key, lines):
    '''Get index of the column specified by header key.

    Inputs:
    key - title of the column
    lines - lines of the input file  

    '''
    header = lines[0].split('\t')
    for i, col in enumerate(header):
        if col == key:
            idx = i
            break
    return idx


def add_uniprot(key, path, params):
    '''Add column of Uniprot ids using the tid column as input.

    Inputs:
    key - title of the column
    path - path to the input file
    params - parameters  

    '''
    infile = open(path, 'r')
    lines = infile.readlines()
    infile.close()
    out = open('_'.join([path,"sed"]) ,'w')
    out.write('%s\tuniprot\n'%lines[0].rstrip('\n') )
    idx = get_idx(key, lines)
    uniprots ={}
    for line in lines[1:]:
        elements = line.split('\t')
        tid = int(elements[idx])
        uniprots[tid] = []
    tidstr = "','".join(map(str, uniprots.keys()))
    print "Looking up uniprot ids for ", len(uniprots.keys()), "proteins."
    data = queryDevice.queryDevice("SELECT distinct tid, accession FROM component_sequences cs JOIN target_components tc ON tc.component_id = cs.component_id WHERE tid IN('%s')"% tidstr, params)
    for tup in data:
        tid = int(tup[0])
        uniprot = tup[1]
        uniprots[tid].append(uniprot)
    uniprots['UNCHECKED'] = 'not assigned'
    for line in lines[1:]:
        elements = line.split('\t')
        tid = int(elements[idx])
        uniprot = uniprots[tid]
        uniprot = ','.join(uniprot)
        out.write("%s\t%s\n"%(line.rstrip('\n'), uniprot ))
    out.close()
    os.system('mv %s %s'% ('_'.join([path,"sed"]), path))


def add_smiles(key, path, params):
    '''Add column of Uniprot ids using the tid column as input.

    Inputs:
    key - title of the column
    path - path to the input file
    params - parameters  

    '''
    infile = open(path, 'r')
    lines = infile.readlines()
    infile.close()
    out = open('_'.join([path,"sed"]) ,'w')
    out.write('%s\tsmiles\n'%lines[0].rstrip('\n') )
    idx = get_idx(key, lines)
    mols ={}
    for line in lines[1:]:
        elements = line.split('\t')
        mol_chembl_id = elements[idx]
        mols[mol_chembl_id] = None
    tidstr = "','".join(map(str, mols.keys()))
    print "Looking up smiles for ", len(mols.keys()), "compounds."
    data = queryDevice.queryDevice("""
		SELECT distinct md.chembl_id, cs.canonical_smiles 
		FROM molecule_dictionary md 
		JOIN compound_structures cs 
		  ON md.molregno = cs.molregno 
		WHERE md.chembl_id IN('%s')"""% tidstr, params)
    for tup in data:
        mol_chembl_id = tup[0]
        smiles = tup[1]
        mols[mol_chembl_id] = smiles
    for line in lines[1:]:
        elements = line.split('\t')
        mol_chembl_id = elements[idx]
        smiles = mols[mol_chembl_id]
        out.write("%s\t%s\n"%(line.rstrip('\n'), smiles ))
    out.close()
    os.system('mv %s %s'% ('_'.join([path,"sed"]), path))



def add_assigned_target(lkp, key, path, params):
    '''Add column of Uniprot ids using the tid column as input.

    Inputs:
    lkp - output of the catch_query function
    key - title of the column
    path - path to infile
    params = parameters from loca.yaml

    '''
    infile = open(path, 'r')
    lines = infile.readlines()
    infile.close()
    out = open('_'.join([path,"sed"]) ,'w')
    out.write('%s\tassigned_target\n'%lines[0].rstrip('\n') )
    idx = get_idx(key, lines)
    for line in lines[1:]:
        elements = line.split('\t')
        chembl_id = elements[idx]
        target = lkp[chembl_id]
        out.write("%s\t%s\n"%(line.rstrip('\n'), target ))
    out.close()
    os.system('mv %s %s'% ('_'.join([path,"sed"]), path))




