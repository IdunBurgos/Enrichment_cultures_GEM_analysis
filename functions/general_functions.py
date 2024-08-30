import pandas as pd
import os

def merge_phyla(phyla):
    if isinstance(phyla,float):
        return phyla
    
    
    phyla_list = phyla.split("_")

    if len(phyla_list)>1:
        return phyla_list[0]
    else:
        return phyla

def read_allmags_data(): 
    """
    Reads and processes MAGs (Metagenome-Assembled Genomes) data from an Excel file.

    This function performs the following steps:
    1. Reads the data from the "Coverage" sheet of an Excel file.
    2. Sets the "MAG" column as the index.
    3. Processes the "Phylum" column to merge names using the `merge_phyla` function.
    4. Performs basic assertions to ensure data integrity.

    Returns:
        pd.DataFrame: A DataFrame containing the processed MAGs data with the "Phylum" column processed in merge_phyla.
    
    Raises:
        AssertionError: If the number of rows is not 101 or if the returned object is not a DataFrame.
    """  
    directory = os 
    directory = os.fsencode("../output/GEMs/GEMs_no_constraints/")

    all_mags_paper = pd.read_excel("../input/files_from_fairdomhub/All_Mags_for_paper_analysis.xlsx",sheet_name="Coverage_vol2")
    all_mags_paper.set_index("MAG",inplace=True)
    all_mags_paper.drop("Coverage (%)",axis=1,inplace=True)

    #mags_fasta = set([str(string)[2:-5] for string in os.listdir(directory)])
    #all_mags_paper = all_mags_paper[all_mags_paper.index.isin(mags_fasta)].copy()
    all_mags_paper = all_mags_paper[all_mags_paper.new_coverage>1].copy()
    
    phylum_series = all_mags_paper["Phylum"].copy()
    all_mags_paper.drop("Phylum",axis=1)
    all_mags_paper["Phylum"] = phylum_series.map(merge_phyla)
    
    assert len(all_mags_paper.index)==73,"There should be 73 rows in this dataframe"
        
    assert isinstance(all_mags_paper,pd.DataFrame), "Something went wrong with processing of dataframe"
    return all_mags_paper


def chebi_selected(interesting_super_classes=['carbohydrate derivatives', 'amino acids and derivatives',
       'nucleotides and derivatives', 'oligosaccharides', 'oligopeptides',
       'simple sugars', 'urea and urea derivatives', 'gases',
       'B-vitamins', 'cofactors', 'fatty acids',
       'carboxylic acids and anions', 'alcohols and aldehydes',
       'lipids'], rxn_based=True):
    
    met_chebi_class = pd.read_csv("../output/met_chebi_class.tsv",index_col=0,sep="\t")
    chebi_interesting = met_chebi_class.copy()
    
    if rxn_based:
        chebi_interesting["reaction"] = chebi_interesting.index.map(lambda x: "R_EX_"+x.split("M_",1)[1]).copy()
        chebi_interesting.set_index("reaction",inplace=True)

    # Get rows from met_chebi_class that match the superclass
    chebi_interesting = chebi_interesting[chebi_interesting["self defined super class"].isin(interesting_super_classes)].copy()
    
    return chebi_interesting

    