"""
    Function: getPfam

    Gets Pfam domains from a list of ChEMBL ids.
    -------------------

    Felix Kruger
    momo.sander@googlemail.com
"""

def getPfam(release):
    import queryDevice
    import os
    import yaml
    import numpy as np

    # Read config file.
    paramFile = open('gla.yaml')
    params = yaml.safe_load(paramFile)
    user = params['user']
    pword = params['pword']
    host = params['host']
    port = params['port']

    # Read infile.
    infile = "data/drugs.tab" 
    infile = open(infile, 'r')
    lines = infile.readlines()
    infile.close()
    chembl_ids = {}
    name_lkp = {}
    for line in lines[1:]:
        elements = line.split('\t')
        chembl_id = elements[1]
        name = elements[4].rstrip()
        name_lkp[chembl_id] = name
        conc = elements[2]
        conc_nm = float(conc) * 1000
        chembl_ids[chembl_id] = conc_nm

    # Query ChEMBL.
    id_string = "','".join(chembl_ids.keys())
    query = """SELECT md.chembl_id, act.molregno, act.activity_id, act.standard_value, act.standard_type, act.assay_id, 
            ass.src_id, 
            ass.tid, ass.relationship_type, 
            td.target_type, td.pref_name
            FROM molecule_dictionary md
            JOIN activities act
            ON md.molregno = act.molregno
            JOIN assays ass
            ON act.assay_id = ass.assay_id
            JOIN target_dictionary td
            ON ass.tid = td.tid
            WHERE md.chembl_id IN('%s')
            AND standard_units = 'nM'
            AND act.potential_duplicate IS NULL
            AND act.standard_relation IN('=', '<')
            AND ass.src_id = 1
            AND standard_type IN('EC50', 'IC50', 'Kd', 'Ki', 'Potency')
            AND relationship_type IN('D', 'H')""" % id_string
    data = queryDevice.queryDevice(query ,release, user, pword, host ,port)

    # Catch query in lkp dictionary.
    lkp = {}
    for ent in data:
        chembl_id = ent[0]
        target = ent[10]
        target_type = ent[9]
        pot = ent[3]
        assay_id =  ent[5]
        if pot > chembl_ids[chembl_id]:
            print 'chucking out:', pot, chembl_id
            continue
        pot = -np.log10(pot)+9
        try:
            lkp[chembl_id][target].append(pot)
        except KeyError:
            try:
                lkp[chembl_id][target] = [pot]
            except KeyError:
                lkp[chembl_id] = {}
                lkp[chembl_id][target] = [pot]

    # Write outfiles. Uses the lkp and name_lkp.
    out = open('data/summary.tab', 'w')
    out_long = open('data/query_results.tab', 'w')
    out_long.write("chembl_id\ttarget\ttarget_type\tconc\tassay_id\tmedian\tn\tname\n")
    out.write('chembl_id\ttarget\tn\tmedian\tn\n')
    for chembl_id in lkp.keys(): 
        for target in lkp[chembl_id].keys():
            n = len(lkp[chembl_id][target])
            median = np.median(lkp[chembl_id][target])
            for pot in lkp[chembl_id][target]:
                out_long.write("%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\n" % (chembl_id, target, target_type, pot, assay_id, median, n, name_lkp[chembl_id] ))
            out.write("%s\t%s\t%s\t%.2f\n" % (chembl_id, target, n ,median)) 
    out.close()
    out_long.close()


if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2: # the program name and the two arguments 
        sys.exit("Must specify path to file, release, user, pword, host, port")
    release = sys.argv[1]
    redThread(release)

