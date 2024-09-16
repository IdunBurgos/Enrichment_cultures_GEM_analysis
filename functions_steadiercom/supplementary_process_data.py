import reframed
import pandas as pd
import matplotlib as mpl

import sys
sys.path.append("../functions/")
sys.path.append("../functions_steadiercom/")
#sys.path.append("../functions_data_analysis/")


import general_functions as general_func
import steadiercom_samples_preprocessing as steadiercom_pre
import colors_community
import colors_MAGs
"""
change_name_dict = {"Isobutyrate":"Isobutyric acid",
"2 methylpropanoic acid":"Isobutyric acid",
"Propionate (n-C3:0)":"Propionic acid",                                                                            
"Butyrate (n-C4:0)":"Butyric acid",
"L-Lactate":"L-Lactic acid",
"Formate":"Formic acid",
"Acetate":"Acetic acid",
"Hexanoate (n-C6:0)":"Caproic acid",
"CO2 CO2":"CO2",
"L-Malate":"L-Malic acid",
"D-Aspartate":"D-Aspartic acid",
"L-Aspartate":"L-Aspartic acid",
"Succinate":"Succinic acid",
"Thymine C5H6N2O2":"Thymine",
"Glycolate C2H3O3":"Glycolic acid",
"R Acetoin C4H8O2":"R Acetoin",
"Protoheme C34H30FeN4O4": "Protoheme",
"N-Acetyl-D-glucosamine(anhydrous)N-Acetylmuramic acid":"N-Acetylmuramic acid",
"R R 2 3 Butanediol C4H10O2":"(R,R)-butane-2,3-diol",
"L Ornithine C5H13N2O2":"L-Ornithine",
"Methanethiol CH4S":"Methanethiol",
"Pentanoate":"Valeric acid",
"Isovalerate":"Isovaleric acid",
"Lactose C12H22O11":"Lactose",
"Urea CH4N2O":"Urea",
"Riboflavin C17H20N4O6":"Riboflavin",
"2 methylbutyraldehyde C5H10O":"2-Methylbutyraldehyde",
"Diacetyl C4H6O2":"Diacetyl"
}

model_uni = reframed.load_cbmodel("../input/universe_bacteria.xml")


# Create new Genus and Family name

all_mags_paper = general_func.read_allmags_data()
#dont_know = all_mags_paper[all_mags_paper["new_coverage"]>10][["Source","Substrate","Family","Genus","new_coverage"]].sort_values(["Source","Substrate"])

all_mags_paper_reduced,total_members_genus = general_func.all_mags_paper_genus(all_mags_paper,prefix=True)
all_mags_paper_reduced,total_members_family = general_func.all_mags_paper_family(all_mags_paper_reduced,prefix=True,combine=False)
genus_groups,mag2genus_dict = general_func.mag2genus(all_mags_paper_reduced)
family_groups,mag2family_dict = general_func.mag2family(all_mags_paper_reduced)




def genus_donor(row):
    if row.donor=="environment":
        return "environment"
    else:
        return mag2genus_dict[row.donor]
    
def genus_receiver(row):
    if row.receiver=="environment":
        return "environment"
    else:
        return mag2genus_dict[row.receiver]
    
def family_donor(row):
    if row.donor=="environment":
        return "environment"
    else:
        return mag2family_dict[row.donor]

def family_receiver(row):
    if row.receiver=="environment":
        return "environment"
    else:
        return mag2family_dict[row.receiver]
    

def met2metname(met):
    met_name = model_uni.metabolites[met].name
    return met_name


def change_name(x):
    if x in change_name_dict.keys():
        return change_name_dict[x]
    else:
        return x

"""

#phyla_lut, unique_phyla, phylum_colors = colors_MAGs.phylum_colors_func()
phyla_lut, unique_phyla, phylum_colors = colors_MAGs.phylum_colors_func()

def process_data(steadier_sample,all_mags_paper_reduced):
    genus_groups,mag2genus_dict = general_func.mag2genus(all_mags_paper_reduced)
    family_groups,mag2family_dict = general_func.mag2family(all_mags_paper_reduced)
    
    
    # Read and process data
    steadier_sample["mass_rate*frequency"] = steadier_sample["mass_rate"]*steadier_sample["frequency"]
    steadier_sample.loc[:,"super_class"] = steadier_sample.compound.map(steadiercom_pre.assign_super_class)

    steadier_sample.compound = steadier_sample.compound.map(steadiercom_pre.met2metname)
    steadier_sample.compound = steadier_sample.compound.map(steadiercom_pre.change_name)
    steadier_sample.compound = steadier_sample.compound.map(lambda x: x.replace("C4H10O2",""))

    # Define groups of different genera and families

    MAGs_steady_com = list(set(list(steadier_sample.donor.values)+list(steadier_sample.receiver.values)))
    MAGs_steady_com.remove("environment")


    steadier_sample.loc[:,"genus_donor"] = steadier_sample.apply(steadiercom_pre.genus_donor,mag2genus_dict=mag2genus_dict,axis=1).copy()
    steadier_sample.loc[:,"genus_receiver"] = steadier_sample.apply(steadiercom_pre.genus_receiver,mag2genus_dict=mag2genus_dict,axis=1).copy()
    steadier_sample.loc[:,"family_donor"] = steadier_sample.apply(steadiercom_pre.family_donor,mag2family_dict=mag2family_dict,axis=1).copy()
    steadier_sample.loc[:,"family_receiver"] = steadier_sample.apply(steadiercom_pre.family_receiver,mag2family_dict=mag2family_dict,axis=1).copy()
    return steadier_sample


def phylum_colors_spec(steadier_sample,all_mags_paper_reduced):
    
    relevant_phyla = all_mags_paper_reduced.Phylum.unique()
    
    genera = list(set(list(steadier_sample.genus_donor.unique()) + list(steadier_sample.genus_receiver.unique())))
    genera.remove("environment")

    phyla = [all_mags_paper_reduced.loc[all_mags_paper_reduced.loc[all_mags_paper_reduced.Genus==genus].index]["Phylum"][0] for genus in genera]
    colors_list = [phyla_lut[phylum] for phylum in phyla]
    
    color_df = pd.DataFrame({"Genus_names":genera,"Phylum":colors_list})
    color_df.set_index("Genus_names",inplace=True)
    return color_df,relevant_phyla,phyla_lut

def color_df_auxotrophies(CD_X_MAGs,all_mags_paper_reduced):
    genus_groups,mag2genus_dict = general_func.mag2genus(all_mags_paper_reduced)
    color_dict = {}
    phyla_specific = []
    for MAG in all_mags_paper_reduced.index:

        phylum = all_mags_paper_reduced.loc[MAG,"Phylum"]
        color_dict[mag2genus_dict[MAG]]=phyla_lut[phylum]
        if MAG in CD_X_MAGs:
            phyla_specific.append(phylum)

    phyla_specific = list(set(phyla_specific))
    colrs_df = pd.Series(color_dict)
    colrs_df.name="Phylum"
    return colrs_df,phyla_specific,phyla_lut
    