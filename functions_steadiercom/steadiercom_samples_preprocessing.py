import reframed
import pandas as pd
import matplotlib as mpl

import sys
sys.path.append("../functions/")
#sys.path.append("../functions_data_analysis/")
import general_functions as general_func


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

## Superclass
met2superclass_dict = pd.read_csv("../output/met_chebi_class.tsv",sep="\t",index_col=0)["self defined super class"].to_dict()

def assign_super_class(x):
    super_class = "other"
    if x in met2superclass_dict.keys():
        super_class = met2superclass_dict[x]
    return super_class    

def genus_donor(row,mag2genus_dict):
    if row.donor=="environment":
        return "environment"
    else:
        return mag2genus_dict[row.donor]
    
def genus_receiver(row,mag2genus_dict):
    if row.receiver=="environment":
        return "environment"
    else:
        return mag2genus_dict[row.receiver]
    
def family_donor(row,mag2family_dict):
    if row.donor=="environment":
        return "environment"
    else:
        return mag2family_dict[row.donor]

def family_receiver(row,mag2family_dict):
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
    
    
