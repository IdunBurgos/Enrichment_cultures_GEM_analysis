#!/usr/bin/env python
# coding: utf-8

# ## Compounds produced by commmunity

print("**COMPOUNDS PRODUCED EXPERIMENTALLY BY COMMUNITY**")
compounds_dict, source_dict,substrate_dict, gas_sheet_dict, community_dict = translation_dicts.translation_dicts()

compounds_dict['Isovalericacid(mg/L)'] = '3mb'
compounds_dict['Hydrogen'] = 'h2'
compounds_dict['CO2'] = 'co2'

compounds_dict_list = {name:[id_] for name,id_ in compounds_dict.items()}
compounds_dict_list["Isobutyricacid(mg/L)"] = ['isobuta','ibt','2mpa']
compounds_dict_list["Isovalericacid(mg/L)"] = ['3mb','3mb','ival']
compounds_dict_list['Caproicacid(mg/L)'] = ['caproic','hxa']
compounds_dict_list["Levulinicacid"] = ['4opa','4oxptn']
compounds_dict_list['Hydrogen'] = ['h2']
compounds_dict_list['CO2'] = ['co2']


# ### Make overview of all compounds measured the enrichment cultures
# 
# The strategy for defining if the compound is produced or not is:
# - Assuming that positive concertation==the compound is produced by the community

community_production_names = {}

# For each source 
for name, id_source in source_dict.items():
    
    # Read and process excel sheet
    ## Get the Excel sheet specifically for the source
    data_df = pd.read_excel("../input/files_from_fairdomhub/enrichment_cultures _data.xlsx",sheet_name=name)
    data_df.dropna(how='all', axis='columns',inplace=True)
    
    data_df = data_df.iloc[:,1:-1]
    data_df.set_index("Sample name ",inplace=True)
    data_df.columns = [col.replace(" ","") for col in data_df.columns]
    
    # For each substrate
    for name_sub,id_sub in substrate_dict.items():
        community_id = id_source+"_"+id_sub
                
        if community_id=="M_A": # Marshland on Avicel is not a combination in our experiments.
            continue
        
        community_production_names[community_id]={}
        data_df_sub_max = data_df[data_df.index.str.contains(community_id)].max()
        
        # For each compound in our analysis
        for compound in compounds_dict.keys():
            if compound not in data_df_sub_max.index:
                continue
            
            # If max concentration is higher than 0 -> assume that it is produced by a community member
            if data_df_sub_max[compound]>0:
                community_production_names[community_id][compound]=1
                
        # These gasses where observed for all enrichment cultures. (for gasses CH4 is not included because it is not a part of the bacterial universal model)
        community_production_names[community_id]["Hydrogen"]=1
        community_production_names[community_id]["CO2"]=1

# ### Save data 

print("\t Save data")

with open("../output/community_production_names.json", "w") as outfile: 
    json.dump(community_production_names, outfile)


with open("../output/compounds_dict_list.json", "w") as outfile: 
    json.dump(compounds_dict_list, outfile)


with open("../output/compounds_dict.json", "w") as outfile: 
    json.dump(compounds_dict, outfile)





