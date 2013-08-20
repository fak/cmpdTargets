import queryDevice
import yaml
import addProperties
import helper
import os
def compound_profile():
    """
    Function: compound_profile

    Retrieves all activities and targets for a list of compounds specified 
    in a tab-separated file.
    -------------------

    Felix Kruger
    momo.sander@googlemail.com
    """
    # Read config file.
    paramFile = open('local.yaml')
    params = yaml.safe_load(paramFile)
    paramFile.close()
    # Read source table.
    (conc_thresholds, assigned_targets, names) =  helper.read_cmpd_table(params['query_name'])
    # Query ChEMBL.
    id_string = "','".join(conc_thresholds.keys())
    query = """SELECT md.chembl_id, 
            act.standard_value, act.assay_id, act.standard_relation, 
            ass.tid, ass.assay_type, 
            td.target_type, td.pref_name,   
	    dcs.pubmed_id
            FROM molecule_dictionary md
            JOIN activities act
            ON md.molregno = act.molregno
            JOIN assays ass
            ON act.assay_id = ass.assay_id
            JOIN target_dictionary td
            ON ass.tid = td.tid
            JOIN docs dcs
            ON ass.doc_id = dcs.doc_id
            WHERE md.chembl_id IN('%s')
            AND standard_units = 'nM'
            AND act.potential_duplicate IS NULL
            AND act.standard_relation IN('=', '<', '>')
            AND ass.src_id = %s
            AND standard_type IN('EC50', 'IC50', 'Kd', 'Ki', 'Potency')
            AND relationship_type IN('D', 'H')""" % (id_string, params['src_id'])
    data = queryDevice.queryDevice(query ,params)
    lkp = helper.capture(data, conc_thresholds)
    helper.write_out(params['query_name'], lkp, names, assigned_targets)
    # Add uniprot ids and canonical smiles.
    path = 'data/%s_results.tab' % params['query_name']
    addProperties.add_uniprot('tid', path, params)
    addProperties.add_smiles('chembl_id', path, params)
    # Quote table columns.
    helper.reformat(path, path + '.tmp')
    os.system("mv %s %s" % (path + '.tmp', path))


if __name__ == '__main__':
    import sys
    if len(sys.argv) != 1: # the program name and the two arguments 
        sys.exit("All arguments are specified in local.yaml")
 
    compound_profile()

