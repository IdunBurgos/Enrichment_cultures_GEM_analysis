#!/usr/bin/env python
# coding: utf-8


# NB: This function does NOT confirm that the metabolite is produced, it just checks if it has the exchange reaction.
def best_candidate(SC_rxns,MAGs,diff_df_select):
    """
    Input:
    SC_rxns - list of exchanged metabolites of community (exchange reactions)
    MAGs - list of MAGs in this community
    diff_df_select - dataframe with increase in rxns and symmetric difference
    """
    rxn_candidates ={}
    
    # For each exchanged compound find the candidate with lowest score in symmetric_diff
    for rxn in SC_rxns:
        rxn_candidates[rxn]={"MAG":None,"symmetric_diff":1000} # Default
        
        # For each community member
        for MAG in MAGs:
            
            # If reaction is not in model it is not a candidate
            if rxn not in GEMs_dict[soft_constraint_selected][MAG].reactions:
                continue
                
            symmetric_diff = diff_df_select.loc[MAG,"symmetric_diff"]  
            
            # If the symmetric difference is smaller than for the previous candidate
            if symmetric_diff<rxn_candidates[rxn]["symmetric_diff"]:
                rxn_candidates[rxn]["MAG"]=MAG
                rxn_candidates[rxn]["symmetric_diff"]=symmetric_diff
                
    return rxn_candidates

print("**SELECT MODELS**")

print("\t Load data...")
compounds_dict, source_dict,substrate_dict, gas_sheet_dict, community_dict = translation_dicts.translation_dicts()
relevant_MAGs_reduced = copy.copy(relevant_MAGs)
relevant_MAGs_reduced.remove('CH15-bin.15') # Had issues with this model for reconstruction with constraints
relevant_MAGs_reduced.remove('CH14-bin.0')

#Load models made without soft constraints


print("\t Load models...")
directory = os.fsencode("../output/GEMs/GEMs_no_constraints/")

GEMs_dict = {"no_constr":{}}
for file in os.listdir(directory):
    
    filename = os.fsdecode(file)
    if filename.endswith(".xml"): 
        GEMs_dict["no_constr"][filename[:-4]]= reframed.load_cbmodel("../output/GEMs/GEMs_no_constraints/"+filename)


#Load models made with soft constraints
directories = {
               "constr0_1":"../output/GEMs/GEMs_intermediate/GEMs_soft_constraints_score_0.1/"
              }

for id_,directory_str in directories.items():
    directory = os.fsencode("../output/GEMs/GEMs_intermediate/GEMs_soft_constraints_score_0.1/")

    GEMs_dict[id_] = {}
    for file in os.listdir(directory):

        filename = os.fsdecode(file)
        if filename.endswith(".xml"): 
            GEMs_dict[id_][filename[:-4]]= reframed.load_cbmodel(directory_str+filename)


#### Find difference between models made with different reconstruction parameters
#- new_rxns_count = number of new reactions
#- symmetric_diff = number of lost reactions + number of new reactions 

print("\t Find difference between models...")
difference_dict ={}
    
for constr_status in directories.keys():

    difference_dict[constr_status] = {}

    difference_dict[constr_status]["new_rxns_count"] = {}
    difference_dict[constr_status]["symmetric_diff"] = {}

    for MAG in relevant_MAGs_reduced:
        model_const = GEMs_dict[constr_status][MAG]
        model_no_constr = GEMs_dict["no_constr"][MAG]

        # number of new reactions in new model
        difference_dict[constr_status]["new_rxns_count"][MAG]=len(set(model_const.reactions)-set(model_no_constr.reactions))
        # reaction symmetric difference
        difference_dict[constr_status]["symmetric_diff"][MAG]=len(set(model_const.reactions).symmetric_difference(set(model_no_constr.reactions)))

diff_df = pd.DataFrame.from_dict({(constr_status,diff_type):difference_dict[constr_status][diff_type] 
                        for constr_status in difference_dict.keys() 
                        for diff_type in difference_dict[constr_status].keys()})

#### Select models: The best candidates to fulfill the community metabolic phenotype
#  Here we focus on the most promising sets of reconstructed models. For this we chose:
#- constr0_1: low soft constraint score will likely lead to addition of fewer reactions
#- no_constr: no soft constraints included. This is the ideal, but not all metabolites were produced.
# Selection criteria: The models that required the least amount of changes (symmetric difference) to acquire the desired phenotype.

# change this if you want to study something different.
soft_constraint_selected = "constr0_1"

#Find best candidates - Based on the fact that we expect certain compounds produced, which models created with soft constraints are the best candidates? Through this code we select just one model for each community and compound produced.
print("\t Select models...")
diff_df_select = diff_df.xs(soft_constraint_selected,axis=1)

# A new key,value pair for the genome-scale metabolic models from 'soft_constraint_selected' or 'no_costr'
GEMs_dict["adapt"]={}

# An overview of the origin of the models in GEMs_dict["adapt"]={}
GEMs_adapt = {}

rxn_candidates_all = {}

# For each community 

for community_name, community_id in community_dict.items():
    # Find community members
    MAGs = MAG2community_id[MAG2community_id[1]==community_id].index.values
    MAGs = [MAG for MAG in MAGs if MAG in relevant_MAGs_reduced] # Only look at 99% best ones.
    
    # Find the exchange reactions representing produced compounds by this community
    soft_constraints = pd.read_csv(f"../output/soft_constraints/SC_{community_id}.tsv", header=None, sep="\t")
    SC_rxns = soft_constraints[0].values
    
    # Find best candidate in community for producing each compound
    rxn_candidates = best_candidate(SC_rxns,MAGs,diff_df_select) 
    rxn_candidates_all[community_id]=rxn_candidates
            
    # Create GEMs adapt
    MAGs_best_candidates = set([candidate_dict["MAG"] for rxn,candidate_dict in rxn_candidates.items()])
    
    for MAG in MAGs:
        
        if MAG in MAGs_best_candidates:
            GEMs_dict["adapt"][MAG]=GEMs_dict[soft_constraint_selected][MAG].copy()
            GEMs_adapt[MAG] = soft_constraint_selected
        else:
            GEMs_dict["adapt"][MAG]=GEMs_dict["no_constr"][MAG].copy()
            GEMs_adapt[MAG] = "no_constr"


# Some candidates are quite different from the original model. Here the symmetric diff shows how many reactions are different between the community members. 

best_candidates = pd.DataFrame.from_dict({(community_id,rxn):rxn_candidates_all[community_id][rxn]  
                                          for community_id in rxn_candidates_all.keys() 
                                          for rxn in rxn_candidates_all[community_id].keys()}).transpose()

##### Save data

GEMs_adapt["CH15-bin.15"]="no_constr"

with open("../output/GEMs/GEMs_intermediate/GEMs_adapt/GEMs_adapt.json", "w") as outfile: 
    json.dump(GEMs_adapt, outfile)

##### Save models

print("\t Save models...")

GEMs_dict["adapt"]["CH15-bin.15"] = GEMs_dict["no_constr"]["CH15-bin.15"].copy()

for MAG,model in GEMs_dict["adapt"].items():
    reframed.save_cbmodel(model,"../output/GEMs/GEMs_intermediate/GEMs_adapt/"+MAG+".xml")