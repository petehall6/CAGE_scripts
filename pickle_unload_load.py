import pickle
import pprint


obj = pickle.load(open('params','rb'))

with open("out_check.txt","a") as f:
    pprint.pprint(obj, stream=f)



""" 


args  = {'GENBANK_LENGTH_EXTENSION': '5000',
 'customer': 'YANG',
 'experiment': 'TEST',
 'gene': 'CRLF2',
 'project_objective': 'CODON CHANGE',
 'species': 'human'}

dict(args)

PROJECT_PARAMS = 'params'

with open(PROJECT_PARAMS, "wb") as fh:
            params = args.copy()
            pickle.dump(params, fh, -1) """