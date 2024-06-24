import numpy as np
from matplotlib.patches import Rectangle
from pycirclize import Circos
from matplotlib.patches import Patch

import matplotlib.pyplot as plt
import reframed
import pandas as pd


from steadiercom_data_processing import mmol2g,find_members,data_ReceiverOrDonor,data_community_abundance_func,data_for_links,data_for_links2
import colors_MAGs as color_func
import general_functions as general_func



## Universal data
all_mags_paper = general_func.read_allmags_data()
chebi_lut, chebi_interesting, chebi_colors_ser = color_func.chebi_rxn_color_func(rxn_based=False)
chebi_lut["other"] = "#eaeaea"
phyla_lut, unique_phyla, phylum_colors = color_func.phylum_colors_func()   

met2superclass_dict = chebi_interesting["self defined super class"].to_dict()

def assign_super_class(row):
    met = row.compound
    
    super_class = "other"
    if met in met2superclass_dict.keys():
        super_class = met2superclass_dict[met]
    return super_class

def data_uptake_prod(steadiercom_crossfeeding_all,abundance,community_id=False,compound_type=False,fluxorrate="flux",only_media=False):

    ## Process input data (add superclass and convert to mass_rate)
    
    if community_id:
        steadiercom_crossfeeding_all_copy = steadiercom_crossfeeding_all.xs(community_id).copy()
        sim_abundance = abundance.xs(community_id).copy()
        
    else:
        steadiercom_crossfeeding_all_copy = steadiercom_crossfeeding_all.copy()
        sim_abundance = abundance.copy()
        
    phyla_groups,mag2phyla_dict,MAGs_steadycom_dict = find_members(steadiercom_crossfeeding_all_copy,all_mags_paper) # NEW
        
    if fluxorrate=="flux":
        steadiercom_crossfeeding_all_copy.loc[:,"mass_rate"] = steadiercom_crossfeeding_all_copy.apply(mmol2g,fluxorrate=fluxorrate,axis=1)

    steadiercom_crossfeeding_all_copy.loc[:,"super_class"]= steadiercom_crossfeeding_all_copy.apply(assign_super_class,axis=1).copy()
    
    
    if compound_type:
        steadiercom_crossfeeding_all_copy.loc[:,"mass_rate"] = steadiercom_crossfeeding_all_copy.apply(lambda x: 0 if x.super_class not in compound_type else x.mass_rate,axis=1)

    
    # DATA CONSUMERS

    data_receiver_df = data_ReceiverOrDonor(steadiercom_crossfeeding_all_copy,"receiver",phyla_groups,mag2phyla_dict,MAGs_steadycom_dict,only_media=only_media)

    # DATA PRODUCERS
    data_donor_df = data_ReceiverOrDonor(steadiercom_crossfeeding_all_copy,"donor",phyla_groups,mag2phyla_dict,MAGs_steadycom_dict,only_media=only_media)
    
    # DATA FOR ABUNDANCE "PLOT
    community_abundance = data_community_abundance_func(all_mags_paper,sim_abundance,community_id,phyla_groups,mag2phyla_dict,MAGs_steadycom_dict)

    
    return data_receiver_df,data_donor_df,community_abundance




def plot_uptake_prod(data_receiver_df,data_donor_df,community_abundance,title=None):
    from matplotlib.lines import Line2D

    fig, axs = plt.subplots(4,gridspec_kw={'height_ratios': [1,3,0.5, 3]},figsize=(20,15))

    # Relative abundance
    community_abundance = community_abundance.reset_index()
    community_abundance["index"]=community_abundance["index"].map(str)
    community_abundance.plot(ax=axs[0],kind="scatter",x="index",y="exp_relative_abundance",ylabel="",xlabel="",color="red")
    #community_abundance.plot(ax=axs[0],kind="scatter", x="index",y="sim_relative_abundance",color="red",label="simulated",style="x",s=50,xlabel="",ylabel="")
    axs[0].set_xticks([])
    axs[0].set_title("Relative abundance",x=-.1,y=0.5,fontsize=15)


    # Producing bars
    data_donor_df = data_donor_df.loc[data_donor_df.sum(axis=1)>1e-20,:].copy()
    data_donor_df.sort_index(axis=1).sort_index(ascending=True).transpose().plot(kind="bar",ax=axs[1],color = chebi_lut,legend=None,stacked=True)
    axs[1].set_xticks([])
    axs[1].set_title("Compound produced \n [g/(h*g$_{bio}$)]",x=-.1,y=0.5,fontsize=15)

    # Phyla colors
    data_ex = data_receiver_df.sort_index(axis=1).transpose()
    xticks_labels = data_ex.index
    xticks = np.arange(0,len(xticks_labels)+1)

    axs[2].set_yticks([0,1])


    y_min, y_max = axs[3].get_ylim()

    axs[2].set_xticks(xticks)
    for i,x in enumerate(xticks_labels):
        patch=Rectangle((i, y_min), 1, y_max-y_min, color=phyla_lut[x[0]])
        axs[2].add_patch(patch)

    axs[2].set_xticks([])
    axs[2].set_yticks([])
    axs[2].set_title("Phylum",x=-.07,y=0.55,fontsize=15)

    plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.1, 
                        hspace=0.05)
    # Consuming bars
    data_receiver_df = data_receiver_df.loc[data_receiver_df.sum(axis=1)>1e-20,:].copy()
    data_receiver_df.sort_index(axis=1).sort_index(ascending=False).transpose().apply(lambda x: (-1)*x).plot(kind="bar",ax=axs[3],color=chebi_lut,legend=None,stacked=True)
    axs[3].set_xticks([])
    axs[3].set_title("Compound consumed \n [g/(h*g$_{bio}$)]",x=-.1,y=0.5,fontsize=15)


    custom_lines = [Patch(color=chebi_lut[super_class], label=super_class) for super_class in  pd.concat([data_donor_df,data_receiver_df]).index.unique()]

    plt.legend(custom_lines, pd.concat([data_donor_df,data_receiver_df]).index.unique(),fontsize=15,loc= (1.01,1.5))
    fig.suptitle(title,fontsize=30)
    plt.show()
    return fig



def circos_plot_process_data(steadiercom_crossfeeding_all,community_id=False,min_flux=0,scale_by_flux=False,compound_type=False,fluxorrate="flux"):
    met2superclass_dict = chebi_interesting["self defined super class"].to_dict()
    
    
    if community_id:
        steadiercom_crossfeeding_all_copy = steadiercom_crossfeeding_all.xs(community_id).copy()
    else:
        steadiercom_crossfeeding_all_copy = steadiercom_crossfeeding_all.copy()
    
    
    phyla_groups,mag2phyla_dict,MAGs_steadycom_dict = find_members(steadiercom_crossfeeding_all_copy,all_mags_paper)
    
    
    if fluxorrate=="flux":
        steadiercom_crossfeeding_all_copy.loc[:,"mass_rate"] = steadiercom_crossfeeding_all_copy.apply(mmol2g,fluxorrate=fluxorrate,axis=1)
    else:
        assert "mass_rate" in steadiercom_crossfeeding_all_copy.columns, "mass_rate not in columns of steadiercom_crossfeeding_all_copy"
    
    steadiercom_crossfeeding = steadiercom_crossfeeding_all_copy.dropna().copy()
    
    
    ## Production data (add superclass and convert to mass_rate)

    if compound_type:
        mets_interesting = chebi_interesting[chebi_interesting["self defined super class"].isin(compound_type)].index.values
        steadiercom_crossfeeding = steadiercom_crossfeeding[steadiercom_crossfeeding.compound.isin(mets_interesting)].copy()

    links,member_index = data_for_links(steadiercom_crossfeeding,phyla_groups,mag2phyla_dict,MAGs_steadycom_dict,min_flux=min_flux,scale_by_flux=scale_by_flux,fluxorrate=fluxorrate)
    #return links,phyla_groups
    

    steadiercom_crossfeeding_all_copy = steadiercom_crossfeeding_all_copy[steadiercom_crossfeeding_all_copy.compound.isin(met2superclass_dict.keys())].copy() # NB here it will exclude compounds that are not found in our classification system. 
    steadiercom_crossfeeding_all_copy.loc[:,"super_class"]= steadiercom_crossfeeding_all_copy.apply(lambda x:met2superclass_dict[x.compound],axis=1).copy()

    # Remove data for non-relevant compound types.
    if compound_type:
        steadiercom_crossfeeding_all_copy.loc[:,"mass_rate"] = steadiercom_crossfeeding_all_copy.apply(lambda x: 0 if x.super_class not in compound_type else x.mass_rate,axis=1)

    data_df = data_ReceiverOrDonor(steadiercom_crossfeeding_all_copy,"donor",phyla_groups,mag2phyla_dict,MAGs_steadycom_dict)
    
    return links,data_df,phyla_groups



def circos_plot_process_data2(steadiercom_crossfeeding_all,community_id=False,min_flux=0,compound_type=False,fluxorrate="flux"):
    met2superclass_dict = chebi_interesting["self defined super class"].to_dict()
    
    
    if community_id:
        steadiercom_crossfeeding_all_copy = steadiercom_crossfeeding_all.xs(community_id).copy()
    else:
        steadiercom_crossfeeding_all_copy = steadiercom_crossfeeding_all.copy()
    
    
    phyla_groups,mag2phyla_dict,MAGs_steadycom_dict = find_members(steadiercom_crossfeeding_all_copy,all_mags_paper)
    
    
    if fluxorrate=="flux":
        steadiercom_crossfeeding_all_copy.loc[:,"mass_rate"] = steadiercom_crossfeeding_all_copy.apply(mmol2g,fluxorrate=fluxorrate,axis=1)
    else:
        assert "mass_rate" in steadiercom_crossfeeding_all_copy.columns, "mass_rate not in columns of steadiercom_crossfeeding_all_copy"
    
    steadiercom_crossfeeding = steadiercom_crossfeeding_all_copy.dropna().copy()
    
    
    ## Production data (add superclass and convert to mass_rate)

    if compound_type:
        mets_interesting = chebi_interesting[chebi_interesting["self defined super class"].isin(compound_type)].index.values
        steadiercom_crossfeeding = steadiercom_crossfeeding[steadiercom_crossfeeding.compound.isin(mets_interesting)].copy()

    links,member_index = data_for_links2(steadiercom_crossfeeding,phyla_groups,mag2phyla_dict,MAGs_steadycom_dict,min_flux=min_flux,fluxorrate=fluxorrate)
    #return links,phyla_groups
    

    return links,phyla_groups





def plot_circos_plot(links,data_df,phyla_groups,color_by_source=False,fontsize=15,title=None):
    ### Plot
    color_dict = chebi_lut
    
    sectors = {phylum:len(MAGs) for phylum,MAGs in phyla_groups.items()}
    name2color = phyla_lut

    circos = Circos(sectors,space=5)
    
    vmax = data_df.sum().max()

    for sector in circos.sectors:
        
        # Add taxonomic group
        track = sector.add_track((60, 65))
        track.axis(fc=name2color[sector.name])
        
        
        # Plot bar track
        bar_track = sector.add_track((65, 100), r_pad_ratio=0.1)
        bar_track.axis()


        data_df_sub = data_df.xs(sector.name,axis=1).transpose().copy()
        data_df_sub.sort_index(inplace=True)
        
        data_df_sub.index = [str(ent) for ent in data_df_sub.index]
        sb_table = bar_track.stacked_bar(
            data_df_sub,
            width=0.3,
            vmax=vmax,
            cmap=color_dict,
            bar_kws=dict(ec="black", lw=0.2),
            label_pos="bottom",
            label_kws=dict(size=0, orientation="horizontal")
        )


    colors_classes = []
    
    for link in links:
        if color_by_source:
            color = name2color[link[0][0]]
            class_ = link[0][0]
        else:
            if link[3] in chebi_colors_ser.keys():
                color = chebi_colors_ser[link[3]]
                class_ = chebi_interesting.loc[link[3],"self defined super class"]
            else:
                color = "#eaeaea"
                class_ = "other" 

        circos.link(link[0],link[1],direction=1,color=color,ec=color,lw=1,arrow_length_ratio=0.15,allow_twist=False)#,alpha=alpha)

        colors_classes.append((color,class_))

    
     ### Plot

    
    fig = circos.plotfig(figsize=(15,15))

    extra_colors = [(chebi_lut[class_],class_) for class_ in data_df[data_df.sum(axis=1)>0].index]
    if color_by_source==False:
        colors_classes.extend(extra_colors)
    
    colors_classes = list(set(colors_classes))
    
    line_handles = [Patch(color=color_code[0], label=color_code[1]) for color_code in colors_classes]
    line_legend = circos.ax.legend(
        handles=line_handles,
        bbox_to_anchor=(0.93, 1.02),
        fontsize=fontsize,
        title="COMPOUNDS",
        handlelength=2,
    )
    plt.title(title,fontsize=30)
    
    return fig

def plot_circos_plot2(links,phyla_groups,color_by_source=False,fontsize=15,title=None):
    ### Plot
    color_dict = chebi_lut

    sectors = {phylum:len(MAGs) for phylum,MAGs in phyla_groups.items()}
    name2color = phyla_lut

    circos = Circos(sectors,space=5)

    for sector in circos.sectors:

        # Add taxonomic group
        track = sector.add_track((95, 100))
        track.axis(fc=name2color[sector.name])

        track = sector.add_track((90, 95))
        track.grid(y_grid_num=None,x_grid_interval=1,alpha=1.0,color="black")

    colors_classes = []

    for link in links:
        if color_by_source:
            color = name2color[link[0][0]]
            class_ = link[0][0]
        
        else:
            if link[3] in chebi_colors_ser.keys():
                color = chebi_colors_ser[link[3]]
                class_ = chebi_interesting.loc[link[3],"self defined super class"]
            else:
                color = "#eaeaea"
                class_ = "other"

        #circos.link(link[0],link[1],direction=1,color=color,ec=color,lw=1)#,alpha=alpha)
        circos.link_line(link[0],link[1],direction=1,color=color,lw=link[2]*200,arrow_height=5,arrow_width = 7,alpha=0.7)


        colors_classes.append((color,class_))
     ### Plot


    fig = circos.plotfig(figsize=(15,15))

    colors_classes = list(set(colors_classes))
    
    line_handles = [Patch(color=color_code[0], label=color_code[1]) for color_code in colors_classes]
    line_legend = circos.ax.legend(
        handles=line_handles,
        bbox_to_anchor=(0.93, 1.02),
        fontsize=fontsize,
        title="COMPOUNDS",
        handlelength=2,
    )
    plt.title(title,fontsize=30)
    
    return fig