"""Script:  loader.py

Load the mapping of Pfam-A domains. The main function is loader, defined at the bottom of the document. Specify release and version of the mapping on command line. eg.: $> python loader.py chembl_15 0_1

Note on variable names: lkp is used to represent dictionaries I was too lazy to find a proper name for.

--------------------
Author:
Felix Kruger
fkrueger@ebi.ac.uk

"""
import csv
import numpy as np



def reformat(infile, outfile):
    '''Quote table columns.

    Inputs:
    infile -- filepath to source table
    outfile -- filepath to output table
    '''
    ifl = open(infile, 'r')
    reader = csv.reader(ifl, delimiter = '\t')
    
    ofl =  open(outfile, 'w')
    writer = csv.writer(ofl, delimiter='\t', quoting=csv.QUOTE_ALL)
    writer.writerows(reader)
    ifl.close()
    ofl.close()


def read_cmpd_table(query_name):
    '''Read table that specifies compounds for which profiles are to be 
       generated.

    Inputs:
    query_name - name of the query and source file 
    '''
    # Read infile.
    infile = "data/%s.tab" % query_name
    infile = open(infile, 'r')
    lines = infile.readlines()
    infile.close()
    conc_thresholds = {}
    names = {}
    assigned_targets = {}
    for line in lines[1:]:
        elements = line.split('\t')
        assigned_target = elements[0]
        chembl_id = elements[1]
        name = elements[4].rstrip()
        names[chembl_id] = name
        conc = elements[2]
        conc_nm = float(conc) * 1000
        conc_thresholds[chembl_id] = conc_nm
        assigned_targets[chembl_id] = assigned_target
    return(conc_thresholds, assigned_targets, names)




def capture(data, conc_thresholds):
    '''Capture query data in a look-up table. 

    Inputs:
    data - output from the SQL query
    conc_thresholds - threshold associated with each compound 

    '''
    lkp = {}
    for ent in data:
        mol_chembl_id = ent[0]
        pot = ent[1]
        assay_id = ent[2]
        std_rel = ent[3]
        tid = ent[4]
        assay_type = ent[5]
        target_type = ent[6]
        pref_name = ent[7]
        pubmed = ent[8]
        if pot > conc_thresholds[mol_chembl_id] or std_rel == '>':
            pot = 4
            #print rel_type, pot
        else:
           pot = -np.log10(pot)+9
        try:
            lkp[mol_chembl_id][pref_name].append((pot, assay_id, target_type, tid, pubmed, assay_type))
        except KeyError:
            try:
                lkp[mol_chembl_id][pref_name] = [(pot, assay_id, target_type, tid, pubmed, assay_type)]
            except KeyError:
                lkp[mol_chembl_id] = {}
                lkp[mol_chembl_id][pref_name] = [(pot, assay_id, target_type, tid ,pubmed, assay_type)]
    return lkp





def write_out(query_name, lkp, names, assigned_targets):
    '''Write out the query results.

    Inputs:
    query_name - name of the query and source file 
    lkp - output of the capture function
    names - dict[chembl_id] = name, output of the read_cmpd_table function
    assigned_targets - dict[chembl_id] = target, output of read_cmpd_table  

    '''
    # Write outfiles. Uses the lkp and name_lkp.
    out = open('data/%s_results.tab' % query_name,  'w')
    out.write("chembl_id\tpref_name\tconc\tmedian\tn\tassay_id\tassay_type\ttype\ttid\tpubmed\tassigned_target\tname\n")
    for chembl_id in lkp.keys():
        for pref_name in lkp[chembl_id].keys():
            n = len(lkp[chembl_id][pref_name])
            median = np.median([x[0] for x in lkp[chembl_id][pref_name]])
            target_type = lkp[chembl_id][pref_name][0][2]
            assigned_target = assigned_targets[chembl_id]
            name = names[chembl_id]
            for pot in lkp[chembl_id][pref_name]:
                pot_v = pot[0]
                assay_id = pot[1]
                pubmed = pot[4]
                tid = pot[3]
                assay_type = pot[5]
                out.write("%(chembl_id)s\t%(pref_name)s\t%(pot_v).2f\t%(median).2f\t%(n)i\t%(assay_id)s\t%(assay_type)s\t%(target_type)s\t%(tid)s\t%(pubmed)s\t%(assigned_target)s\t%(name)s\n" % locals())
    out.close()



