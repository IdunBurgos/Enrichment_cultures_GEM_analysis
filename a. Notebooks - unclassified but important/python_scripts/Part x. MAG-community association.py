#!/usr/bin/env python
# coding: utf-8

# # Data showing MAG-community association



print("**CREATE MAG-COMMUNITY ASSOCIATION**")

compounds_dict, source_dict,substrate_dict, gas_sheet_dict, community_dict = translation_dicts.translation_dicts()

# ### Read data and models
# **Read data**
all_Mags_for_paper_analysis = general_func.read_allmags_data()

MAG2community_id = pd.Series({MAG:community_dict[community_id] for community_id,MAGs in all_Mags_for_paper_analysis.groupby(["Source","Substrate"]).groups.items() for MAG in MAGs})
MAG2community_id


# ## Save data
print("\t Save data")
MAG2community_id.to_csv("../output/MAG2community_id.tsv",sep="\t",header=False)




