#!/usr/bin/env python
# coding: utf-8

print("**MAKING SOFT CONSTRAINTS**")

print("\t read data")
compounds_dict, source_dict,substrate_dict, gas_sheet_dict, community_dict = translation_dicts.translation_dicts()
all_Mags_for_paper_analysis = general_func.read_allmags_data()
MAG2community_id = pd.read_csv("../output/MAG2community_id.tsv",sep="\t",header=None,index_col=0)

with open("../output/relevant_MAGs_99.txt") as text_file:
    relevant_MAGs = text_file.read().split("\n")
relevant_MAGs = [string.replace("\t","") for string in relevant_MAGs]

with open("../output/community_production_names.json") as text_file:
    community_production_names = json.load(text_file)

with open("../output/compounds_dict_list.json") as text_file: 
    compounds_dict_list = json.load(text_file)

with open("../output/compounds_dict.json") as text_file: 
    compounds_dict = json.load(text_file)
    
print("\t Load models..")
directory = os.fsencode("../output/GEMs/GEMs_no_constraints/")
GEMs_dict = {}
for file in os.listdir(directory):
    
    filename = os.fsdecode(file)
    if filename.endswith(".xml"): 
        GEMs_dict[filename[:-4]]= reframed.load_cbmodel("../output/GEMs/GEMs_no_constraints/"+filename)
        
        
##### Defining soft constraints

#- Excluding elements that can already be produced by a community member in complete media
#- Focus on 99% most abundant species

##### Find the producers (from the 99% most abundant species)

print("\t determine soft constraints")
## Find the community they belong to
MAG2community_id_most_abundant = MAG2community_id[MAG2community_id.index.isin(relevant_MAGs)]
enrich_groups_top = MAG2community_id_most_abundant.groupby(1).groups# top 99

## Find producers in COMPLETE media and growth=0

# For each community set the default production to False
comm_producers_top = {community_id:{MAG:{compound:False for compound in compounds.keys()} for MAG in enrich_groups_top[community_id]} for community_id,compounds in community_production_names.items()}
for community_id,MAGs in enrich_groups_top.items():
    for MAG in MAGs:
        if MAG not in GEMs_dict.keys():
            continue
        model = GEMs_dict[MAG].copy()
        
        complete_env = reframed.Environment.from_model(model)      
        
        for compound in community_production_names[community_id].keys():
            mets = compounds_dict_list[compound]

            ex_rxns = ["R_EX_"+met+"_e" for met in mets if "R_EX_"+met+"_e" in model.get_exchange_reactions()]
            
            if len(ex_rxns)==0:
                continue
            
            sol = reframed.FVA(model,constraints=complete_env,reactions=ex_rxns)
            total_prod = sum([value[1] for value in sol.values()])
            comm_producers_top[community_id][MAG][compound] = total_prod>1e-6
            
dfs_community_count_top = {community_id:pd.DataFrame(producers_dict).sum(axis=1).rename("top99%") for community_id,producers_dict in comm_producers_top.items()}

##### Define soft constraints

soft_constraints_new = {}

for community_id,producers_df in dfs_community_count_top.items():
    # Find compounds not produced
    not_produced = producers_df[producers_df==0].index
    
    soft_constraints_new[community_id]={}
    
    for compound in not_produced:
        soft_constraints_new[community_id]["R_EX_"+compounds_dict[compound]+"_e"]=1
        
##### Media for gapfilling during reconstruction
#We know that some of these compounds are produced when a specific electron donor is present. 

media_db = pd.read_csv("https://raw.githubusercontent.com/cdanielmachado/carveme/master/carveme/data/input/media_db.tsv",sep="\t")

lb_02 = media_db[media_db.medium=="LB[-O2]"].copy()
lb_02.medium = lb_02.medium.map(lambda x:x.replace("LB[-O2]","LB_extend"))
lb_02.loc[-1] = ["LB_extend", "Additional elements for MCCA production", "lac__L","L-Lactate"] 
lb_02.loc[-1] = ["LB_extend", "Additional elements for MCCA production", "etoh","Ethanol"] 
lb_02.loc[-1] = ["LB_extend", "Additional elements for MCCA production", "ac","Acetate"] 
lb_02.loc[-1] = ["LB_extend", "Additional elements for MCCA production", "ppa","Propionate (n-C3:0)"] 
lb_02.loc[-1] = ["LB_extend", "Additional elements for MCCA production", "xyl__D","D-Xylose"] 
lb_02.loc[-1] = ["LB_extend", "Additional elements for MCCA production", "lcts","Lactose"] 

lb_02.reset_index(drop=True,inplace=True)


##### Save data

print("\t Save data")
for community_id, dict_ in soft_constraints_new.items():
    pd.DataFrame(pd.Series(dict_)).to_csv("../output/soft_constraints/SC_"+community_id+".tsv",
                                          sep="\t",
                                          header=False,
                                          index_label=False)
lb_02.to_csv("../output/soft_constraints/SC_media_db.tsv",
                                         sep="\t",
                                         index=False,
                                         index_label=False)

##### Check that everything is as it should be
print("\t Assertions")

for community_id, dict_ in soft_constraints_new.items():
    SC_media_test_old = pd.read_csv("../output_30_08_24/soft_constraints/SC_"+community_id+".tsv",
                                          sep="\t",
                                          header=None)
    SC_media_test_new = pd.DataFrame(pd.Series(dict_)).reset_index()
    SC_media_test_new.columns=[0,1]
    assert  SC_media_test_old[SC_media_test_old==SC_media_test_new].shape==SC_media_test_old.shape
    assert  SC_media_test_old[SC_media_test_old==SC_media_test_new].shape==SC_media_test_new.shape
    
old_lb_02 = pd.read_csv("../output_30_08_24/soft_constraints/SC_media_db.tsv",
                                         sep="\t")

assert old_lb_02.shape==old_lb_02[old_lb_02==lb_02].shape
assert lb_02.shape==old_lb_02[old_lb_02==lb_02].shape