import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import colormaps
import matplotlib.pyplot as plt

from pycirclize import Circos
import reframed
import pandas as pd


from steadiercom_samples_processing import preprocessing_func, data_ReceiverOrDonor,data_community_abundance_func,data_for_links


import sys
sys.path.append('../functions')
import colors_MAGs as color_func
import general_functions as general_func


chebi_lut, chebi_interesting, chebi_colors_ser = color_func.chebi_rxn_color_func(rxn_based=False)
chebi_lut["other"] = "#eaeaea"
phyla_lut, unique_phyla, phylum_colors = color_func.phylum_colors_func() 

all_mags_paper = general_func.read_allmags_data()
cmap = colormaps.get("gist_ncar")

def circos_plot_process_data(steadiercom_samples,all_mags_paper =all_mags_paper,community_id=False,min_flux=0,compound_type=False,min_frequency=0):
    
    # Add super_class and mass_rate*frequency
    steadiercom_samples_preprocessed = preprocessing_func(steadiercom_samples)
    

    # Select data for community
    if community_id:
        steadiercom_samples_preprocessed_copy = steadiercom_samples_preprocessed[steadiercom_samples_preprocessed.community==community_id].copy()
    else:
        steadiercom_samples_preprocessed_copy = steadiercom_samples_preprocessed.copy()
    
    # find members
    steadiercom_crossfeeding_only_producers = steadiercom_samples_preprocessed_copy[(steadiercom_samples_preprocessed_copy.donor!="environment")]
    steadiercom_crossfeeding_only_consumers = steadiercom_samples_preprocessed_copy[(steadiercom_samples_preprocessed_copy.receiver!="environment")]
    members =  set(list(steadiercom_crossfeeding_only_producers.donor.unique())+list(steadiercom_crossfeeding_only_consumers.receiver.unique()))
    
    # sort the members according to coverage
    members = all_mags_paper[all_mags_paper.index.isin(members)].sort_values("new_coverage",ascending=False).index
    
    # remove all data related to environment
    steadiercom_samples_preprocessed_copy = steadiercom_samples_preprocessed_copy[(steadiercom_samples_preprocessed_copy.donor!="environment") & (steadiercom_samples_preprocessed_copy.receiver!="environment")]
    
    # filter out lower frequencies 
    steadiercom_samples_preprocessed_copy = steadiercom_samples_preprocessed_copy[steadiercom_samples_preprocessed_copy.frequency>min_frequency].copy()
    
    ## Production data (add superclass and convert to mass_rate)
    if compound_type:
        mets_interesting = chebi_interesting[chebi_interesting["self defined super class"].isin(compound_type)].index.values
        steadiercom_crossfeeding = steadiercom_samples_preprocessed_copy[steadiercom_samples_preprocessed_copy.compound.isin(mets_interesting)].copy()
    else:
        steadiercom_crossfeeding = steadiercom_samples_preprocessed_copy.copy()
    
    # Sort crossfeeding data according to the 
    steadiercom_crossfeeding.donor = steadiercom_crossfeeding.donor.astype("category")
    steadiercom_crossfeeding.donor = steadiercom_crossfeeding.donor.cat.set_categories(members)
    steadiercom_crossfeeding.receiver = steadiercom_crossfeeding.receiver.astype("category")
    steadiercom_crossfeeding.receiver = steadiercom_crossfeeding.receiver.cat.set_categories(members)
    steadiercom_crossfeeding = steadiercom_crossfeeding.sort_values(["donor","receiver"],ascending=[False,True]).copy()
    

    links = data_for_links(steadiercom_crossfeeding,members,min_flux=min_flux)
    return links,members



def plot_circos_plot(links,members,fontsize=15,title=None):
    ### Plot
    color_dict = chebi_lut

    sectors = {member:1 for member in members}
    #name2color = phyla_lut

    circos = Circos(sectors,space=5)
    
    i = 0
    
    index_color = np.arange(0.1,1.0,0.9/len(circos.sectors))
    
    for sector in circos.sectors:
        
        # Add color for group
        track = sector.add_track((95, 100))
        
        track.axis(fc=cmap(index_color[i]))
        i+=1

        
    colors_classes = []

    for link in links:

        if link[2] in chebi_colors_ser.keys():
            color = chebi_colors_ser[link[2]]
            class_ = chebi_interesting.loc[link[2],"self defined super class"]
        else:
            color = "#eaeaea"
            class_ = "other" 

        circos.link(link[0],link[1],direction=1,color=color,ec=color,arrow_length_ratio=0.15,allow_twist=False)


        colors_classes.append((color,class_))
     ### Plot


    fig = circos.plotfig(figsize=(15,15))

    colors_classes = list(set(colors_classes))
    
    line_handles = [Patch(color=color_code[0], label=color_code[1],alpha=0.8) for color_code in colors_classes]
    line_legend = circos.ax.legend(
        handles=line_handles,
        bbox_to_anchor=(0.93, 1.02),
        fontsize=fontsize,
        title="COMPOUNDS",
        handlelength=2,
    )
    plt.title(title,fontsize=30)
    
    mag2color= dict(zip(members,index_color))
    
    return fig,mag2color


def data_uptake_prod(steadiercom_samples,all_mags_paper=all_mags_paper,community_id=False,compound_type=False,only_media=False,min_frequency=0.0):
    
    
    steadiercom_samples_preprocessed = preprocessing_func(steadiercom_samples)
    
    
    
    # filter out lower frequencies 
    steadiercom_samples_preprocessed = steadiercom_samples_preprocessed[steadiercom_samples_preprocessed.frequency>min_frequency].copy()
    
    
    ## Process input data (add superclass and convert to mass_rate)
    
    if community_id:
        
        steadiercom_samples_preprocessed_copy = steadiercom_samples_preprocessed[steadiercom_samples_preprocessed.community==community_id].copy()
        
    else:
        steadiercom_samples_preprocessed_copy = steadiercom_samples_preprocessed.copy()

        
    if compound_type:
        steadiercom_samples_preprocessed_copy.loc[:,"mass_rate*frequency"] = steadiercom_samples_preprocessed_copy.apply(lambda x: 0 if x.super_class not in compound_type else x["mass_rate*frequency"],axis=1)

    
    # DATA CONSUMERS

    data_receiver_df = data_ReceiverOrDonor(steadiercom_samples_preprocessed_copy,"receiver",only_media=only_media)

    # DATA PRODUCERS
    data_donor_df = data_ReceiverOrDonor(steadiercom_samples_preprocessed_copy,"donor",only_media=only_media)
    
    # DATA FOR ABUNDANCE "PLOT
    community_abundance = data_community_abundance_func(steadiercom_samples_preprocessed_copy,all_mags_paper,community_id=community_id)
    
    
    
    ## Sort dataframes
    
    data_receiver_df = data_receiver_df.loc[:,community_abundance.index].copy()
    data_receiver_df = data_receiver_df.sort_index().copy()
    
    data_donor_df = data_donor_df.loc[:,community_abundance.index].copy()
    data_donor_df = data_donor_df.sort_index().copy()
    

    return data_receiver_df,data_donor_df,community_abundance


def plot_uptake_prod(data_receiver_df,data_donor_df,community_abundance,title=None):


    fig, axs = plt.subplots(4,gridspec_kw={'height_ratios': [1,3,0.5, 3]},figsize=(20,15))

    # Relative abundance
    community_abundance = community_abundance.reset_index()

    community_abundance.plot(ax=axs[0],kind="scatter",x="MAG",y="exp_relative_abundance",ylabel="",xlabel="",color="red")
    axs[0].set_xticks([])
    axs[0].set_title("Relative abundance",x=-.1,y=0.5,fontsize=15)


    # Producing bars
    data_donor_df = data_donor_df.loc[data_donor_df.sum(axis=1)>1e-20,:].copy()
    data_donor_df.transpose().plot(kind="bar",ax=axs[1],color = chebi_lut,legend=None,stacked=True)
    axs[1].set_xticks([])
    axs[1].set_xlabel('')
    axs[1].set_title("Compound produced \n [g/(h*g$_{bio}$)]",x=-.1,y=0.5,fontsize=15)

    # Phyla colors
    xticks_labels = data_receiver_df.columns
    xticks = np.arange(0,len(xticks_labels)+1)
    
    index_color = np.arange(0.1,1.0,0.9/len(xticks_labels))
    
    axs[2].set_yticks([])


    axs[2].set_xticks(xticks)
    for i,x in enumerate(xticks_labels):
        patch=Rectangle((i, 0),1, 1, color=cmap(index_color[i]))
        
        

        axs[2].add_patch(patch)

    axs[2].set_xticks([])
    plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.1, 
                        hspace=0.05)
    axs[2].set_xticks([])
   
    # Consuming bars
    data_receiver_df = data_receiver_df.loc[data_receiver_df.sum(axis=1)>1e-20,:].copy()
    data_receiver_df.transpose().apply(lambda x: (-1)*x).plot(kind="bar",ax=axs[3],color=chebi_lut,legend=None,stacked=True)
    axs[3].set_xticks([])
    axs[3].set_xlabel('')
    axs[3].set_title("Compound consumed \n [g/(h*g$_{bio}$)]",x=-.1,y=0.5,fontsize=15)


    custom_lines = [Patch(color=chebi_lut[super_class], label=super_class,alpha=0.8) for super_class in  pd.concat([data_donor_df,data_receiver_df]).index.unique()]

    plt.legend(custom_lines, pd.concat([data_donor_df,data_receiver_df]).index.unique(),fontsize=15,loc= (1.01,1.5))
    fig.suptitle(title,fontsize=30)
    plt.show()
    
    mag2color= dict(zip(xticks_labels,index_color))
    return fig,mag2color