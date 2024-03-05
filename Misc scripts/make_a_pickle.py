import pickle
import pprint

#*Change the following parameters

def make_a_pickle():

    args  = {'GENBANK_LENGTH_EXTENSION': '5000',
    'customer': 'PMH',
    'experiment': 'TEST1', 
    'gene': 'PETEGENE',
    'project_objective': 'KO',
    'species': 'MEATPOPCICLE'}

    dict(args)

    PROJECT_PARAMS = 'params_test'

    with open(PROJECT_PARAMS, "wb") as fh:
                params = args.copy()
                pickle.dump(params, fh, -1) 

def dump_a_pickle():
           
    pickle_dic = pickle.load(open('params_test', 'rb'))
    
    
    with open("pickle_check.txt", "a") as f_handle:
        pprint.pprint(pickle_dic, stream = f_handle)
        
#make_a_pickle()

dump_a_pickle()