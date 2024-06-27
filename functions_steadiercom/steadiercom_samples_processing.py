from molmass import Formula
import pandas as pd
import reframed
import numpy as np
import sys
sys.path.append('../functions')

import translation_dicts
import colors_MAGs as color_func
import general_functions as general_func


all_mags_paper = general_func.read_allmags_data()  
chebi_lut, chebi_interesting, chebi_colors_ser = color_func.chebi_rxn_color_func(rxn_based=False)
chebi_lut["other"] = "#eaeaea"

met2superclass_dict = chebi_interesting["self defined super class"].to_dict()

def assign_super_class(row):
    met = row.compound
    
    super_class = "other"
    if met in met2superclass_dict.keys():
        super_class = met2superclass_dict[met]
    return super_class

def mag2genus(steadiercom_sample,all_mags_paper):  
   
    # Find members of this dataset 
    MAGs_steady_com = set(list(steadiercom_sample[steadiercom_sample.donor!="environment"].donor.values)+list(steadiercom_sample[steadiercom_sample.receiver!="environment"].receiver.values)) 
    # Find the genus for the members

    genus_groups = all_mags_paper[all_mags_paper.index.isin(MAGs_steady_com)].groupby("Genus").groups
    mag2genus_dict = {mag:genus for genus,mags in genus_groups.items() for mag in mags}
    
    return genus_groups,mag2genus_dict,MAGs_steady_com

def preprocessing_func(steadiercom_samples):
    
    steadiercom_samples_preprocessed = steadiercom_samples.copy()
    
    steadiercom_samples_preprocessed["mass_rate*frequency"] = steadiercom_samples_preprocessed.apply(lambda row: row.mass_rate*row.frequency,axis=1)
    
    steadiercom_samples_preprocessed["super_class"] = steadiercom_samples_preprocessed.apply(assign_super_class,axis=1)
    
    return steadiercom_samples_preprocessed


def data_ReceiverOrDonor(steadiercom_samples_preprocessed,ReceiverOrDonor,only_media=False):
    """
    ReceiverOrDonor (string): Either 'receiver' or 'donor'
    """
    steadiercom_only_producers = steadiercom_samples_preprocessed[(steadiercom_samples_preprocessed.donor!="environment")]
    steadiercom_only_consumers = steadiercom_samples_preprocessed[(steadiercom_samples_preprocessed.receiver!="environment")]
    
    members =  set(list(steadiercom_only_producers.donor.unique())+list(steadiercom_only_consumers.receiver.unique()))
    
    if only_media:
        if ReceiverOrDonor=="receiver":
            steadiercom_crossfeeding_from_media = steadiercom_samples_preprocessed[steadiercom_samples_preprocessed.donor=="environment"].copy()
        elif ReceiverOrDonor=="donor":
            steadiercom_crossfeeding_from_media = steadiercom_samples_preprocessed[steadiercom_samples_preprocessed.receiver=="environment"].copy()
    else:
        steadiercom_crossfeeding_from_media = steadiercom_samples_preprocessed.copy()
    
    
    # For each receiver/super_class couple -> sum across axis
    data_receiverOrdonor_df = steadiercom_crossfeeding_from_media.groupby(["super_class",ReceiverOrDonor]).sum()["mass_rate*frequency"].unstack().fillna(0)
    data_receiverOrdonor_df = data_receiverOrdonor_df.drop("environment",axis=1)
    
    ## Add data for missing members:
    missing_members = set(members)-set(data_receiverOrdonor_df.columns)
    for member in missing_members:
        receiver_or_donor = member
        
        data_receiverOrdonor_df[receiver_or_donor]=np.zeros(len(data_receiverOrdonor_df.index))
    
    return data_receiverOrdonor_df

def data_community_abundance_func(steadiercom_samples_preprocessed,all_mags_paper,community_id=False):

    steadiercom_only_producers = steadiercom_samples_preprocessed[(steadiercom_samples_preprocessed.donor!="environment")]
    steadiercom_only_consumers = steadiercom_samples_preprocessed[(steadiercom_samples_preprocessed.receiver!="environment")]
    
    members =  set(list(steadiercom_only_producers.donor.unique())+list(steadiercom_only_consumers.receiver.unique()))
    
    compounds_dict, source_dict,substrate_dict, gas_sheet_dict, community_dict = translation_dicts.translation_dicts()
    substrate_dict = {value:key for key,value in substrate_dict.items()}
    source_dict = {'M': 'Marshland', 'CD': 'Compost_Digestate', 'CM': 'Cow_Manure'}


    # Find coverage
    if community_id:
        community_coverage = all_mags_paper[(all_mags_paper.Source==source_dict[community_id.split("_")[0]]) & (all_mags_paper.Substrate==substrate_dict[community_id.split("_")[1]])]["new_coverage"].copy()
    else:
        community_coverage = all_mags_paper["new_coverage"].copy()
    
    ## PROCESS EXPERIMENTAL DATA
    # Find experimental relative abundance
    community_abundance = community_coverage.map(lambda x: x/community_coverage.sum())
    community_abundance.name="exp_relative_abundance"
    
    
    # Reduce set to the ones in our dataset
    
    community_abundance = community_abundance.sort_values(ascending=False).copy()
    community_abundance = community_abundance[community_abundance.index.isin(members)]
    community_abundance = pd.DataFrame(community_abundance)

    
   
    return community_abundance

def data_for_links(steadiercom_crossfeeding_only,members,min_flux=0):

    links = []
    
    member_range = {member:{"upper_donor":[],"upper_receiver":[]} for member in members}

    sum_fluxes = steadiercom_crossfeeding_only[steadiercom_crossfeeding_only["mass_rate*frequency"]>min_flux]["mass_rate*frequency"].sum()
    
    for index,row in steadiercom_crossfeeding_only.iterrows():
        link = []

        # Donor
        donor =row.donor

        # Receiver
        receiver =row.receiver

        flux_g = row["mass_rate*frequency"]

        
        if flux_g<min_flux:
            continue

        lower_receiver_i = 0 if len(member_range[receiver]["upper_receiver"]) == 0 else member_range[receiver]["upper_receiver"][-1]
        upper_receiver_i = lower_receiver_i + flux_g/sum_fluxes
        
        lower_donor_i = 0 if len(member_range[donor]["upper_donor"]) == 0 else member_range[donor]["upper_donor"][-1]
        upper_donor_i = lower_donor_i + flux_g/sum_fluxes
        

        member_range[receiver]["upper_receiver"].append(upper_receiver_i)

        member_range[donor]["upper_donor"].append(upper_donor_i)
        

        link.append((receiver,lower_receiver_i,upper_receiver_i))
        link.append((donor,1-lower_donor_i,1-upper_donor_i))        

        link.append(row.compound)
        links.append(link)
        
        
    return links