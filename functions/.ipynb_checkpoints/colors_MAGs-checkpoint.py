import general_functions
import pandas as pd
import seaborn as sns
import matplotlib.colors as mcolors
from matplotlib import colors



def phylum_colors_func(): 
    """
        Returns:
        phyla_lut (dict): A dictionary mapping each unique phylum (str) to a color (RGB tuple).
        unique_phyla (ndarray): An array of unique phylum names.
        phylum_colors (pd.Series): A pandas Series mapping each MAG in the dataset to its corresponding color.
    
    Usage Example:
        phyla_lut, unique_phyla, phylum_colors = phylum_colors_func()
    """
    
    all_mags_paper = general_functions.read_allmags_data()

    all_phyla = all_mags_paper["Phylum"].values
    unique_phyla = all_mags_paper["Phylum"].dropna().unique()
    #phyla_pal= sns.color_palette("tab20",15)
    phyla_pal = [ darken_color("#41CD05", amount=0.1), 
    darken_color("#CC0000", amount=0.07), 
    darken_color("#FB50A0", amount=0.07), 
    darken_color("#00caca",amount=0.15), 
    darken_color("#0000ff",amount=0.05), 
    darken_color("#a000f0", amount=0.02), 
    "#A6BC09", 
    "#122070", 
    darken_color("#ffff00",amount=0.05), 
    "#744700", 
    "#AF1C72", 
    darken_color("#F55A00",amount=0.05),
    "#38761D" 
        ]

    #phyla_pal = colors.ListedColormap( cs)
    phyla_lut = dict(zip(map(str, unique_phyla),phyla_pal))

    all_mags_paper_reduced = all_mags_paper["Phylum"]
    phylum_colors = pd.Series(all_phyla,index=all_mags_paper_reduced.index).map(phyla_lut)
    
    phylum_colors.name = "Phylum"
    return phyla_lut,unique_phyla,phylum_colors

def cazy_colors_func():
    
    """
        Returns:
        cazy_lut (dict): A dictionary mapping the number of polysaccharides (str) to a color (RGB tuple).
        unique_numbers (ndarray): An array of number of polysaccharides an organism is likely to degrade.
        cazy_colors (pd.Series): A pandas Series mapping each MAG in the dataset to its corresponding color.
    
    Usage Example:
        cazy_lut,unique_numbers,cazy_colors = cazy_colors_func()
    """
    
    product = pd.read_csv("../input/files_from_fairdomhub/product.tsv",index_col="genome",sep="\t")
    
    product_count = product.filter(regex="CAZy*").sum(axis=1)
    product_count = product_count.map(str)

    unique_numbers = ['0','1','2','3','4']
    cazy_pal= sns.color_palette("Greens",5)

    cazy_lut = dict(zip(map(str, unique_numbers),cazy_pal))
    cazy_lut['0']=(1.0,1.0,1.0) 

    cazy_colors = pd.Series(product_count,index=product_count.index).map(cazy_lut)
    #cazy_colors = cazy_colors.loc[cazy_colors.index.isin(GEMs_dict.keys())]
    
    cazy_colors.name = "CAZy"
    return cazy_lut,unique_numbers,cazy_colors


def darken_color(hex_color, amount=0.4):

    rgb = mcolors.hex2color(hex_color)

    darkened_rgb = [max(0, c - amount) for c in rgb]

    return mcolors.to_hex(darkened_rgb)


def chebi_rxn_color_func(rxn_based=True,selected_super_classes=False):
    
    """
    Args:
        met_based (bool): If colors should be based on metabolite name or reaction name (default=Fales (reaction))
    
    Returns:
        chebi_lut (dict): A dictionary mapping each interesting superclass (str) to a color (RGB tuple).
        chebi_interesting (pd.DataFrame): The filtered DataFrame containing only the interesting superclasses.
        chebi_colors_ser (pd.Series): A pandas Series mapping each reaction ID to its corresponding color.
    
    Usage Example:
        chebi_lut, chebi_interesting, chebi_colors_ser = chebi_rxn_color_func()
    """

    if rxn_based:
        chebi_interesting = general_functions.chebi_selected()
    else:
        chebi_interesting = general_functions.chebi_selected(rxn_based=False)
        
    interesting_super_classes = chebi_interesting["self defined super class"].unique()
    """
    # OLD COLORS
    chebi_pal=[
        "#FF85A1",  # Vivid Bubblegum Pink
        darken_color("#FFAD60",0.1),  # Vivid Peach
        "#FFE680",  # Vivid Light Yellow
        "#8FCFFF",  # Vivid Light Blue
        "#C2A2FF",  # Vivid Light Purple
        "#9BE7A8",  # Vivid Light Green
        "#FFD19B",  # Vivid Light Orange
        "#BAE1FF",  # Pastel Blue    
        "#FFA3B1",  # Vivid Pink
        "#BAFFC9",  # Pastel Green
        darken_color("#E0BAFF",0.1),  # Pastel Purple
        "#FFD1DC",   # Pastel Pink Blush
        "#FFB3BA",  # Pastel Pink
        "#FFFFBA"  # Pastel Yellow
    ]
    
    # NEWER VERSION
    chebi_pal=[
    "#FF85A1",  # carbohydrates and derivatives
    "#ddaa7b",  # amino acids
    "#fbdc5d",  # nucleosides
    "#a7dfb6",  # oligosaccharides
    
    "#FFD19B",  # oligopeptides -FFD19B
    "#d5f588",  # simple sugars -d5f588
    "#FFFFBA",  # urea 
    "#b1fbf6",  # gasses
    darken_color("#adb0e8",0.1),  # B-vitamins
    "#ffd1fd",   # fatty acids
    "#FFB3BA",  # carboxylic acids
    "#C3BBF3",  # alcohols and aldehydes C3BBF3
    "#8FCFFF" # phospholipids    
    ]
    
    """
    

    # New colors
    chebi_lut = {"carbohydrate derivatives":"#e5fa79",  # carbohydrates and derivatives  
 "amino acids and derivatives":"#ddaa7b",  # amino acids
 "nucleotides and derivatives":"#fbdc5d",  # nucleosides
 "oligosaccharides":"#a7dfb6",  # oligosaccharides
 "oligopeptides":"#FFD19B",  # oligopeptides -FFD19B
 "simple sugars":"#bfffa4",  # simple sugars -d5f588
 "urea and urea derivatives":"#FFFFBA",  # urea 
 "gases":"#b1fbf6",  # gasses
 "B-vitamins":"#adb0e8",  # B-vitamins
 "fatty acids":"#ffd1fd",   # fatty acids
 "carboxylic acids and anions":"#FFB3BA",  # carboxylic acids
 "alcohols and aldehydes":"#C3BBF3",  # alcohols and aldehydes C3BBF3
 "lipids":"#8FCFFF", # phospholipids    
 "cofactors":"#FF85A1" 
}
    
    #chebi_pal = [darken_color(color, amount=0.0) for color in chebi_pal]

    #chebi_lut = dict(zip(map(str, interesting_super_classes),chebi_pal))
    chebi_interesting.loc[:,"colors"]=chebi_interesting["self defined super class"].map(lambda x:chebi_lut[x])
    chebi_colors_ser = pd.Series(chebi_interesting["colors"].values,index=chebi_interesting.index)
    
    if rxn_based:
        chebi_prefix="R_EX_"
    else:
        chebi_prefix="M_"
        
    
    assert chebi_interesting.loc[chebi_prefix+"glc__D_e","self defined super class"]=="simple sugars"
    assert chebi_interesting.loc[chebi_prefix+"ala__D_e","self defined super class"]=="amino acids and derivatives"
    assert chebi_interesting.loc[chebi_prefix+"adn_e","self defined super class"]== 'nucleotides and derivatives'
    
    
    assert chebi_colors_ser.loc[chebi_prefix+"glc__D_e"]==chebi_lut["simple sugars"]
    assert chebi_colors_ser.loc[chebi_prefix+"ala__D_e"]==chebi_lut["amino acids and derivatives"]
    assert chebi_colors_ser.loc[chebi_prefix+"adn_e"]==chebi_lut['nucleotides and derivatives']
    
    if selected_super_classes:
        chebi_lut_reduced = {key:value for key,value in chebi_lut.items() if key in selected_super_classes}
        chebi_interesting_reduced = chebi_interesting[chebi_interesting["self defined super class"].isin(selected_super_classes)].copy()
        chebi_colors_ser_reduced = chebi_colors_ser[chebi_colors_ser.index.isin(chebi_interesting_reduced.index)]
    else:
        chebi_lut_reduced = chebi_lut
        chebi_interesting_reduced = chebi_interesting
        chebi_colors_ser_reduced = chebi_colors_ser
    return chebi_lut_reduced,chebi_interesting_reduced,chebi_colors_ser_reduced







