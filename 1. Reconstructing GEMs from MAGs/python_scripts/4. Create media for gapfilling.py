print("**CREATE MEDIA FOR GAPFILLING**")

compounds_dict, source_dict,substrate_dict, gas_sheet_dict, community_dict = translation_dicts.translation_dicts()
relevant_MAGs.remove("CH14-bin.0")

all_mags_paper = general_func.read_allmags_data()

##### Load models made without soft constraints
print("\t Load models...")

GEMs_dict = GEMs_dict3

##### Define environment
syncon_environments = MAG_environments.community_syncon_environments()

#### Modify Environment to support community growth
##### FBA growth predictions

MAG2community_id.columns=["community_id"]

community_groups = MAG2community_id.groupby(by="community_id").groups

FBA_growth = {}
for community_id, MAGs in community_groups.items(): 
    FBA_growth[community_id]={}
    for MAG in MAGs:
        if MAG in relevant_MAGs:
            model = GEMs_dict[MAG]
            syncon_environments[community_id].apply(model,inplace=True,exclusive=True,warning=False)
            sol = reframed.FBA(model)

            if sol is None:
                FBA_growth[community_id][MAG]=None
            else:
                FBA_growth[community_id][MAG]=sol

# Add the source and substrate to this data
growth_community_df = pd.concat([MAG2community_id,pd.Series({GEM:sol.fobj for community_id,GEM_sol_dict in FBA_growth.items() for GEM,sol in GEM_sol_dict.items()})],axis=1)
# Change from float to False or positive
growth_community_df["Grows"] = growth_community_df[0].map(lambda x:x>1e-6)
# Drop the growth float column
growth_community_df.drop(0,axis=1,inplace=True)

# When considering top 99 of members all but Cow_Manure on xylan have a growing community member

MAG_can_grow = growth_community_df[growth_community_df.Grows].index

##### FVA prediction of bacteria that can survive in the media at obj_frac=0

# 1. Find compounds produced by growing community members
# 2. Filter by CHEBI class

# Run FVA for exchange reactions of growing community members
print("\t Determine community production...")

FVA_production = {}

for community_id, MAGs in community_groups.items(): 
    
    FVA_production[community_id]={}
    
    for MAG in MAGs:
        
        # If MAG is among the ones who cannot grow -> continue
        if MAG not in MAG_can_grow:
            continue
        model = GEMs_dict[MAG]
        
        # Apply medium
        syncon_environments[community_id].apply(model,inplace=True,exclusive=True,warning=False)
        # Find FVA solution and obj_frac=0
        FVA_production[community_id][MAG] = reframed.FVA(model,reactions=model.get_exchange_reactions(),obj_frac=0.0)
        

# Find compounds that are being produced

FVA_production_copy = FVA_production.copy()
FVA_production_copy = {community:{MAG:[rxn for rxn,sol in FVA_production_copy[community][MAG].items() if sol[1]>1e-6] for MAG in FVA_production_copy[community].keys()}
                       for community in FVA_production_copy.keys()}

# Combine the results from each community member into community level

community_prod = {community_name:[] for community_name in FVA_production_copy.keys()}

for community_name,mag_prod in FVA_production_copy.items():
    for MAG,rxns in mag_prod.items():
        community_prod[community_name].extend(rxns)
    
    community_prod[community_name] = list(set(community_prod[community_name]))

#### Filter by chebi_class
print("\t Filter results...")

# Some compounds are not interesting for us when it comes to exchange
ignore_classes = ["other","inorganic ions and atoms","oligopeptide","simple sugars","cellodextrin","carbohydrate derivative","carbohydrate acid","oligosaccharides"]

met_chebi_class = pd.read_csv("../output/met_chebi_class.tsv",sep="\t",index_col=0)

met_chebi_class_reduced = met_chebi_class[~met_chebi_class["self defined super class"].isin(ignore_classes)].copy()


met_chebi_class_reduced.loc["M_glc__D_e"]= met_chebi_class.loc["M_glc__D_e"]
met_chebi_class_reduced.loc["M_xyl__D_e"]= met_chebi_class.loc["M_xyl__D_e"]


met_chebi_class_dict = met_chebi_class_reduced["chebi class"].to_dict()

community_prod_dfs= {community_name: pd.DataFrame({"rxns":[rxn for rxn in rxns if "M_"+rxn[5:] in met_chebi_class_dict.keys() ],
                                                   "chebi_class":[met_chebi_class_dict["M_"+rxn[5:]] for rxn in rxns if "M_"+rxn[5:] in met_chebi_class_dict.keys() ]})
                                                  for community_name,rxns in community_prod.items()}

##### Make new environments

model_uni = reframed.load_cbmodel("../input/universe_bacteria.xml")

substrate_dict = {'X': 'Xylan', 'A': 'Avicel', 'P': 'PASC'}
source_dict = {'M':'Marshland soil','CD':'Compost and Digestate', 'CM':'Cow manure'}

syncon_media = {}
for community_id in MAG2community_id.community_id.unique():
    
    # Mets from original syncon environment
    mets_syncon = [rxn[5:-2] for rxn in syncon_environments[community_id].keys()]
    
    # Mets produced by the community members
    mets_produced = [rxn[5:-2] for rxn in community_prod_dfs[community_id].rxns.values]
    mets_syncon.extend(mets_produced)

    syncon_media[community_id]=set(mets_syncon)

media_dfs = []

for community_id,compounds in syncon_media.items():
    media_id_list = community_id.split("_")

    source = media_id_list[0]
    substrate = media_id_list[1]

    # NB! Some of the compounds they used are not in the universal model!
    compounds = [met for met in compounds if "M_"+met+"_e" in model_uni.metabolites]
    compounds_names = [model_uni.metabolites["M_"+met+"_e"].name for met in compounds]
    
    media_df = pd.DataFrame({"medium":[community_id for met in compounds],
               "description": ["Media + produced compounds in "+source_dict[source]+" on " +substrate_dict[substrate] for met in compounds],
               "compound":compounds,
               "name":compounds_names})
    media_dfs.append(media_df)

print("\t Save data...")
    
media_total_df = pd.concat(media_dfs)
media_total_df.to_csv("../output/gapfill_media/gapfill_media.tsv",sep="\t",index=None)