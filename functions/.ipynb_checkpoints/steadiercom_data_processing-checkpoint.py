from molmass import Formula
import pandas as pd
import reframed
import numpy as np
import translation_dicts

model_uni = reframed.load_cbmodel("../input/universe_bacteria.xml")

def mmol2g(row,fluxorrate="flux"):
    molweight = Formula(model_uni.metabolites[row.compound].metadata["FORMULA"]).mass
    flux_g = molweight*row[fluxorrate]/1000
    return flux_g

def find_members(steadiercom_crossfeeding_all,all_mags_paper):  
   
    # Find members of this dataset 
    MAGs_steady_com = set(list(steadiercom_crossfeeding_all.donor.dropna().values)+list(steadiercom_crossfeeding_all.receiver.dropna().values))
    MAGs_steadycom_dict = {MAG:MAG.split("_")[0]+"-"+MAG.split("_")[1]+"."+MAG.split("_")[2] for MAG in  MAGs_steady_com}
    MAGs_steady_com = [MAG.split("_")[0]+"-"+MAG.split("_")[1]+"."+MAG.split("_")[2] for MAG in  MAGs_steady_com]
     
    # Find the phyla the members of the 

    phyla_groups = all_mags_paper[all_mags_paper.index.isin(MAGs_steady_com)].groupby("Phylum").groups
    mag2phyla_dict = {mag:phyla for phyla,mags in phyla_groups.items() for mag in mags}
    
    return phyla_groups,mag2phyla_dict,MAGs_steadycom_dict


def data_ReceiverOrDonor(steadiercom_crossfeeding_all_copy,ReceiverOrDonor,phyla_groups,mag2phyla_dict,MAGs_steadycom_dict,only_media=False):
    """
    ReceiverOrDonor (string): Either 'receiver' or 'donor'
    """
    
    members = set(list(steadiercom_crossfeeding_all_copy.donor.dropna().unique())+list(steadiercom_crossfeeding_all_copy.receiver.dropna().unique()))
    
    if only_media:
        if ReceiverOrDonor=="receiver":
            steadiercom_crossfeeding_from_media = steadiercom_crossfeeding_all_copy[steadiercom_crossfeeding_all_copy.donor.isna()].copy()
        elif ReceiverOrDonor=="donor":
            steadiercom_crossfeeding_from_media = steadiercom_crossfeeding_all_copy[steadiercom_crossfeeding_all_copy.receiver.isna()].copy()
    else:
        steadiercom_crossfeeding_from_media = steadiercom_crossfeeding_all_copy.copy()
    
    data_receiverOrdonor = {}
    members_covered = []
    
    # For each receiver/super_class couple -> sum across axis
    for (receiver_or_donor_steadycom,super_class),value in steadiercom_crossfeeding_from_media.groupby([ReceiverOrDonor,"super_class"]).sum()["mass_rate"].items():
        members_covered.append(receiver_or_donor_steadycom)
        
        # Get data about receiver
        receiver_or_donor = MAGs_steadycom_dict[receiver_or_donor_steadycom]
        phylum_receiver_or_donor = mag2phyla_dict[receiver_or_donor]
        receiver_or_donor_index = list(phyla_groups[phylum_receiver_or_donor]).index(receiver_or_donor)
        

        # Save data
        if (phylum_receiver_or_donor,receiver_or_donor_index,receiver_or_donor) in data_receiverOrdonor.keys():
            data_receiverOrdonor[(phylum_receiver_or_donor,receiver_or_donor_index,receiver_or_donor)][super_class] = value
        else:
            data_receiverOrdonor[(phylum_receiver_or_donor,receiver_or_donor_index,receiver_or_donor)]={}
            data_receiverOrdonor[(phylum_receiver_or_donor,receiver_or_donor_index,receiver_or_donor)][super_class] = value
    data_receiverOrdonor_df = pd.DataFrame(data_receiverOrdonor).fillna(0)
    
    
    
    ## Add data for missing members:

    missing_members = set(members)-set(members_covered)
    for member in missing_members:
        receiver_or_donor = MAGs_steadycom_dict[member]
        phylum_receiver_or_donor = mag2phyla_dict[receiver_or_donor]
        receiver_or_donor_index = list(phyla_groups[phylum_receiver_or_donor]).index(receiver_or_donor)
        
        data_receiverOrdonor_df[(phylum_receiver_or_donor,receiver_or_donor_index,receiver_or_donor)]=np.zeros(len(data_receiverOrdonor_df.index))
    
    return data_receiverOrdonor_df




def data_community_abundance_func(all_mags_paper,sim_abundance,community_id,phyla_groups,mag2phyla_dict,MAGs_steadycom_dict):

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
    community_abundance = community_abundance[community_abundance.index.isin(mag2phyla_dict.keys())]
    # Change index to same starndards as for the bar data (Firmicutes,1)
    community_abundance = community_abundance.reset_index()
    community_abundance["index"] = community_abundance.MAG.map(lambda x:(mag2phyla_dict[x],list(phyla_groups[mag2phyla_dict[x]]).index(x)))
    community_abundance.set_index("index")
    
    ## PROCESS SIMULATION DATA
    sim_abundance.columns=["sim_relative_abundance"]
    sim_abundance = sim_abundance[sim_abundance.index.isin(MAGs_steadycom_dict.keys())].copy()
    sim_abundance["MAG"]=sim_abundance.index.map(lambda x: MAGs_steadycom_dict[x])
    sim_abundance_dict = sim_abundance.set_index("MAG")["sim_relative_abundance"].to_dict()
    
    # COMBINE DATA
    community_abundance = community_abundance[community_abundance.MAG.isin(sim_abundance_dict.keys())]
    community_abundance["sim_relative_abundance"] = community_abundance.MAG.map(lambda x: sim_abundance_dict[x])
    community_abundance = community_abundance.set_index("index").sort_index()[["exp_relative_abundance","sim_relative_abundance"]]


    return community_abundance



def data_for_links(steadiercom_crossfeeding,phyla_groups,mag2phyla_dict,MAGs_steadycom_dict,min_flux=0,scale_by_flux=0,fluxorrate="flux"):
    
    links = []
    
    member_index = []

    for index,row in steadiercom_crossfeeding.iterrows():
        link = []

        # Donor
        donor = MAGs_steadycom_dict[row.donor]
        phylum_donor = mag2phyla_dict[donor]
        donor_index = list(phyla_groups[phylum_donor]).index(donor)

        member_index.append((phylum_donor,donor_index,donor))
        # Receiver
        receiver = MAGs_steadycom_dict[row.receiver]
        phylum_receiver = mag2phyla_dict[receiver]
        receiver_index = list(phyla_groups[phylum_receiver]).index(receiver)
        member_index.append((phylum_receiver,receiver_index,receiver))
        
        flux_g = row.mass_rate
        #molweight = Formula(model_uni.metabolites[row.compound].metadata["FORMULA"]).mass
        #flux_g = molweight*row[fluxorrate]/1000

        if scale_by_flux:

            if flux_g<min_flux:
                continue

            if flux_g>0.01:
                ratio = 1
            
            else:
                ratio = flux_g/0.01

            lower_donor = donor_index + 0.48 - 0.4*ratio
            upper_donor = donor_index + 1 - 0.48 + 0.4*ratio

            lower_receiver = receiver_index + 0.48 - 0.4*ratio
            upper_receiver = receiver_index + 1 - 0.48 + 0.4*ratio
            
            link.append((phylum_donor,lower_donor,upper_donor))
            link.append((phylum_receiver,lower_receiver,upper_receiver))
            
        else:
            link.append((phylum_donor,donor_index+0.48,donor_index+0.52))
            link.append((phylum_receiver,receiver_index+0.48,receiver_index+0.52))


        link.append(flux_g)
        link.append(row.compound)
        links.append(link)
        
        member_index = list(set(member_index))
        
    return links,member_index



def data_for_links2(steadiercom_crossfeeding,phyla_groups,mag2phyla_dict,MAGs_steadycom_dict,min_flux=0,fluxorrate="flux"):
    
    links = []
    
    member_index = []

    for index,row in steadiercom_crossfeeding.iterrows():
        link = []

        # Donor
        donor = MAGs_steadycom_dict[row.donor]
        phylum_donor = mag2phyla_dict[donor]
        donor_index = list(phyla_groups[phylum_donor]).index(donor)

        member_index.append((phylum_donor,donor_index,donor))
        # Receiver
        receiver = MAGs_steadycom_dict[row.receiver]
        phylum_receiver = mag2phyla_dict[receiver]
        receiver_index = list(phyla_groups[phylum_receiver]).index(receiver)
        member_index.append((phylum_receiver,receiver_index,receiver))
        
        flux_g = row.mass_rate
        #molweight = Formula(model_uni.metabolites[row.compound].metadata["FORMULA"]).mass
        #flux_g = molweight*row[fluxorrate]/1000


        if flux_g<min_flux:
            continue

        if flux_g/0.1<1e-2:
            ratio = 1e-2
        elif flux_g/0.1>0.5:
            ratio = 0.5
        else:
            ratio = flux_g/0.1

        donor_index = donor_index + 0.5
        
        receiver_index = receiver_index + 0.5
       
        link.append((phylum_donor,donor_index))
        link.append((phylum_receiver,receiver_index))

        link.append(ratio)
        link.append(row.compound)
        links.append(link)
        
        member_index = list(set(member_index))
        
    return links,member_index

