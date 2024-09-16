#!/usr/bin/env python
# coding: utf-8

# # Define classes of metabolites




print("**DEFINING METABOLITE CLASSES**")

# ### Load universal model
model_uni = reframed.load_cbmodel("../input/universe_bacteria.xml")

# ### Find the best way to map exchange to classes

# #### Get Chebi ids associated with extracellular metabolites
print("\t Get chebi IDs from model metadata")

chebi_ids = {}
all_mets_chebi =[] 

for met in model_uni.metabolites:

    if model_uni.metabolites[met].compartment=="C_e":

        #ex_met.append(met)
        if "chebi" in model_uni.metabolites[met].metadata.keys():
            #chebi_in_met.append(met)
            all_mets_chebi.append(met)
            chebi_ids[met]=model_uni.metabolites[met].metadata["chebi"]
        else:
            chebi_ids[met]=None



# **Find all chebi ids for the exchange reactions in all models (this also includes the ones that do not have chebi_id)**

super_dict = chebi_ids.copy()
super_dict["M_but_e"]= ['CHEBI:17968']
super_dict["M_ppa_e"] = ['CHEBI:17272']
super_dict["M_14glucan_e"]= ['CHEBI:15444']
super_dict["M_2h3mp_e"] = ['CHEBI:133085']
super_dict["M_2mba_e"] = ['CHEBI:37070']
super_dict["M_2mpa_e"] = ['CHEBI:16135']
super_dict["M_3mba_e"] =['CHEBI:28484']
super_dict["M_4abzglu_e"] =['CHEBI:60903']
super_dict["M_4hpro_DC_e"] =['CHEBI:16231']
super_dict["M_val__D_e"] =['CHEBI:27477']
super_dict["M_tyr__D_e"] = ['CHEBI:28479']
super_dict["M_leu__D_e"] = ['CHEBI:28225']
super_dict["M_isocap_e"] = ['CHEBI:74903'] 

all_mets = [met[2:-2] for met in super_dict.keys()]

# ## Define substrate classes based on ChEBI
# 
# It's necessary to define some classes from CHEBI that we are interested in. 
# 
# Strategy:
# 
# - Define main classes. NB: Some compounds might fit into several classes
#     - Define a hierarchy to avoid placing some compounds in a general class
#     - Find the chebi class by recursively walking through the map.


really_large_classes = collections.OrderedDict({
                ## Really Bigg classes
                "organonitrogen compounds":{
                    "CHEBI:35352":"organonitrogen compound"},
                "organosulfur compounds":{
                    "CHEBI:33261":"organosulfur compound"},
                "organophosphorus compounds":{
                    "CHEBI:25710":"organophosphorus compound"},
                "carboxylic acids and anions":{
                    "CHEBI:29067":"carboxylic acid anion",
                    "CHEBI:33575":"carboxylic acid"},
                "carbohydrate derivatives":{
                    "CHEBI:63299":"carbohydrate derivative"}
    })
large_classes = collections.OrderedDict({
                ## Bigg classes
                "alcohols and aldehydes":{
                    "CHEBI:15734": "primary alcohol",
                    "CHEBI:35681":"secondary alcohol",
                    "CHEBI:17478":"aldehyde"},
                "oligosaccharides":{
                    "CHEBI:50699":"oligosaccharide (undefined)"},
                "cofactors":{
                    "CHEBI:23357":"cofactor",
                    "CHEBI:33892":"iron coordination entity"}

})
medium_classes = collections.OrderedDict({
                "aromatic compounds":{
                    "CHEBI:26195":"polyphenol"
                },
                "alkaloids":{
                    "CHEBI:22315":"alkaloid"
                },
                "polyamine":{
                    "CHEBI:88061":"polyamine"
                },
                "lipids":{
                    "CHEBI:16247":"phospholipid"},
                "steroid":{
                    "CHEBI:35341":"steroid"
                },
                "gases":{
                    "CHEBI:138675":"gas molecular entity"}
                
})

main_classes = collections.OrderedDict({
                "amino acids and derivatives":{
                    "CHEBI:37022":"amino-acid anion",
                    "CHEBI:33709":"amino acid",
                    "CHEBI:83821":"amino acid derivative"
},
                "oligopeptides":{                    
                    "CHEBI:25676":"oligopeptide"},
    
                "fatty acids":{
                    "CHEBI:58954":"straight-chain saturated fatty acid anion",
                    "CHEBI:58956":"branched-chain saturated fatty acid anion"},


                "carboxylic acids and anions":{
                    "CHEBI:33576":"sulfur-containing carboxylic acid"},

                ## Carbohydrates
                "simple sugars":{
                    "CHEBI:35381":"monosaccharide",
                    "CHEBI:36233": "disaccharide"},
                "carbohydrate derivatives":{
                    "CHEBI:23639":"deoxy sugar",
                    "CHEBI:33720":"carbohydrate acid", 
                },
                "oligosaccharides":{
                    "CHEBI:22590":"arabinan",
                    "CHEBI:37163":"glucan (undefined)"},
                "nucleotides and derivatives":{
                    "CHEBI:18282":"nucleobase", 
                    "CHEBI:33838":"nucleoside",
                    "CHEBI:26401":"purines",
                    "CHEBI:39447":"pyrimidines",
                    "CHEBI:25608":"nucleoside phosphate"},

                "B-vitamins":{
                   "CHEBI:75769": "B vitamin"},

                "other vitamins":{
                    "CHEBI:12777":"vitamin A",
                    "CHEBI:33234":"vitamin E",
                    "CHEBI:27300":"vitamin D",
                    "CHEBI:176783":"vitamin C",
                    "CHEBI:28384":"vitamin K"},
                "inorganic ions and atoms":{
                    "CHEBI:24835":"inorganic ion",
                    "CHEBI:25585": "nonmetal atom"},

             
                "urea and urea derivatives":{
                    "CHEBI:47857":"ureas",
                    "CHEBI:16199":"urea"},
    
                "aromatic compounds":{
                    "CHEBI:33853":"phenols",
                    "CHEBI:27338": "xylene",
                    "CHEBI:27024":"toluenes"},
                "other":{
                    "CHEBI:26191":"polyol",
                    "CHEBI:23217": "cholines",
                    "CHEBI:24828": "indoles",
                    "CHEBI:26188":"polyketide"}
               })


# **Recursive function to find first parent matching with main classes**
print("\t Assign compound class")


def find_main_class(chebi_id, main_classes):
    # If chebi_id is the id of a  main class -> return value
    if chebi_id in main_classes:
        return chebi_id  

    entity = ChebiEntity(chebi_id)
    parents = [rel.get_target_chebi_id() for rel in entity.get_outgoings() if rel.get_type() == "is_a"]
    
    # If we have reached the end of the graph
    if len(parents)==0:
        return None
    
    
    for parent in parents:
        result = find_main_class(parent, main_classes)  
        if result is not None:
            return result  

    return None

met_chebi_dict = {} 

for main_class_dict_nested in [main_classes,medium_classes,large_classes,really_large_classes]:
    main_class_dict = collections.OrderedDict()
    for d in main_class_dict_nested.values():
        for k, v in d.items():  
            main_class_dict[k]=v
    
    for met_id,chebi_list in super_dict.items():
        if chebi_list is None:
            continue
        
        if met_id in met_chebi_dict.keys():
            continue

        for chebi_id in chebi_list:
            main_class = find_main_class(chebi_id,main_class_dict.keys())
            if main_class is not None:
                met_chebi_dict[met_id]=main_class
                break


# #### Set chebi ids for new metabolites from new models
iIB746 = reframed.load_cbmodel("../input/curated_models/iIB746.xml") # Make sure that this is up to date!

xylo_gluc = [rxn for rxn in iIB746.get_exchange_reactions() if "Q" in rxn and rxn.split("_")[2].isupper()]
ara_xyl = [rxn for rxn in iIB746.get_exchange_reactions() if "A" in rxn and rxn.split("_")[2].isupper()]

cellulose = [rxn for rxn in iIB746.get_exchange_reactions() if "cell" in rxn]
cellulose.append("R_EX_cell6_e")

xylan = [rxn for rxn in iIB746.get_exchange_reactions() if "xyla" in rxn]

new_classes_rxns = {"CHEBI:18233":xylo_gluc,"CHEBI:28427":ara_xyl,"CHEBI:3523":cellulose,"CHEBI:60938":xylan}

new_classes = {"oligosaccharides":{"CHEBI:60938":"glucuronoxylan","CHEBI:28427":"arabinoxylan","CHEBI:18233":"xyloglucan","CHEBI:3523":"cellodextrin"}}

for chebi_id,rxns in new_classes_rxns.items():
    for rxn in rxns:
        met_id = rxn.replace("R_EX_","M_")
        met_id = met_id.replace("EX_","M_")
        met_chebi_dict[met_id]=chebi_id


# #### Combine data as dataframe
# **First we need a dict to map chebi id to class and the self defined super class**
chebi_classes = collections.defaultdict(set)

for class_dict in [medium_classes,large_classes,really_large_classes,main_classes,new_classes]:
    for met_super_class_name,sub_class_dict in class_dict.items():
        
        for chebi_id, chebi_class_name in sub_class_dict.items():  

            chebi_classes[chebi_id]={"chebi class":chebi_class_name,"self defined super class":met_super_class_name}


# **Make a df for all the metabolites we were able to classify**
met_chebi_df = pd.DataFrame(met_chebi_dict.values(),index=met_chebi_dict.keys())
met_chebi_df.columns=["chebi id"]


# **Combine data**
met_chebi_df["chebi class"]=met_chebi_df["chebi id"].map(lambda x:chebi_classes[x]["chebi class"])
met_chebi_df["self defined super class"]=met_chebi_df["chebi id"].map(lambda x:chebi_classes[x]["self defined super class"])

# **Add remaining compounds into 'other'**
not_in_chebi = [met for met in super_dict.keys() if met not in met_chebi_df.index.values ]

for met in not_in_chebi:
    met_chebi_df.loc[met] = [None, "not classified", "other"]


# **Fix some gas molecules that were assigned the wrong class** 
met_chebi_df.loc["M_h2s_e"] = ["CHEBI:138675","gas molecular entity","gases"]
met_chebi_df.loc["M_h2_e"] = ["CHEBI:138675","gas molecular entity","gases"]
met_chebi_df.loc["M_no2_e"] = ["CHEBI:138675","gas molecular entity","gases"]
met_chebi_df.loc["M_no_e"] = ["CHEBI:138675","gas molecular entity","gases"]
met_chebi_df.loc["M_so3_e"] = ["CHEBI:138675","gas molecular entity","gases"]
met_chebi_df.loc["M_o2_e"] = ["CHEBI:138675","gas molecular entity","gases"]
met_chebi_df.loc["M_o2s_e"] = ["CHEBI:138675","gas molecular entity","gases"]
met_chebi_df.loc["M_n2_e"] = ["CHEBI:138675","gas molecular entity","gases"]
met_chebi_df[met_chebi_df["self defined super class"]=="gases"]


#### Save data
print("\t Save data")
met_chebi_df.to_csv("../output/met_chebi_class.tsv",sep="\t")




