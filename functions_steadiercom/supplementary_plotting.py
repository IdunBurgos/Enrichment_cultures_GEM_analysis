import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.gridspec as gridspec
import seaborn as sns


communityid2name_dict = {"CD_A":"Compost digestate - Avicel",
"CD_P":"Compost digestate - PASC",
"CD_X":"Compost digestate - Xylan",
"CM_A":"Cow manure - Avicel",
"CM_P":"Cow manure - PASC",
"CM_X":"Cow manure - Xylan",
"M_P":"Marshland soil - PASC",
"M_X":"Marshland soil - Xylan",
}




def plot_interactions(steadier_sample,cmap_lut,selected_compounds=False,interesting_superclasses=False,fig_title=False,min_number=0,ylim_top=100,inset_plot=True): #,count_or_sum="count"
    
    steadier_sample_pro = steadier_sample.copy()
    
    if interesting_superclasses:
        steadier_sample_pro = steadier_sample_pro[steadier_sample_pro.super_class.isin(interesting_superclasses)].copy()
    
    
    cmap_local = colors.ListedColormap([cmap_lut[community] for community in sorted(steadier_sample_pro.community.unique())])    
    steady_com_count = steadier_sample_pro.groupby(["community","compound"]).sum()["mass_rate*frequency"]
    

    # Only allow members of selected categories to be involved (avoid minerals etc)
    if selected_compounds:
        steady_com_count = steady_com_count.reset_index()[steady_com_count.reset_index().compound.isin(chebi_interesting.index)].set_index(["community","compound"])
    steady_com_count.index.names = ("Community","Compound")
    steady_com_count = pd.DataFrame(steady_com_count)

    # prepare dataframe for plot
    all_compounds = steady_com_count.reset_index().groupby("Compound").sum()["mass_rate*frequency"]
    all_compounds=all_compounds[all_compounds>min_number].index  

    steady_com_count = steady_com_count.loc[(slice(None),list(all_compounds)),:]
    steady_com_count_for_plot = steady_com_count["mass_rate*frequency"].unstack(level=0).fillna(0)
    
    #Sort by the sum of the columns
    steady_com_count_for_plot = steady_com_count_for_plot.transpose()[steady_com_count_for_plot.sum(axis=1).sort_values().index]\
                                                                .transpose()
    steady_com_count_for_plot = steady_com_count_for_plot.transpose().reset_index()
    steady_com_count_for_plot["Community"]=steady_com_count_for_plot["Community"].map(lambda x:communityid2name_dict[x])
    steady_com_count_for_plot = steady_com_count_for_plot.set_index("Community").transpose()
    
    
    
    if inset_plot:  
        ax = steady_com_count_for_plot.plot(kind="bar",figsize=(10,5),cmap=cmap_local,fontsize=12, rot=90,legend=False)
    
        plt.xlabel("")
    
        plt.savefig("Figures/interactions_bar_temporary.png",bbox_inches='tight',dpi=300)

        plt.close()


        fig, ax = plt.subplots()
        im = plt.imread("Figures/interactions_bar_temporary.png")

        steady_com_count_for_plot.plot(ax=ax,kind="bar",figsize=(20,10),cmap=cmap_local,fontsize=12, rot=90,ylim=[0,ylim_top])
    else:
        ax = steady_com_count_for_plot.plot(kind="bar",figsize=(20,10),cmap=cmap_local,fontsize=12, rot=90,legend=False)
    
        
    plt.legend(fontsize=12)
    plt.xlabel('Compounds', fontsize=15)

    plt.ylabel("Interaction rate \n [mg/(h*g$_{bio}$)]",fontsize=15)
    
    if inset_plot:
        newax = fig.add_axes([0.15, 0.3, 0.7, 0.55], zorder=10)
        newax.imshow(im)
        newax.axis('off')


    plt.savefig("Figures/interactions_"+str(fig_title)+".png",bbox_inches='tight',dpi=500)
    plt.show()
    
    
    return steady_com_count_for_plot



def make_receiver_frequency_plots(ax1,ax2,interaction_strength_specific,interaction_strength_all,class_name,get_data=False):
    palette={"all compounds":"grey",
             "carboxylic acids":"red",
             "amino acids, nucleotides, and derivatives":"tab:blue",
             "alcohols and aldehydes":"tab:green"}

    ## Scatter plot
    ax1 =interaction_strength_all.plot(ax=ax1,y="frequency",x="relative_abundance",kind="scatter",label="all compounds",color="grey")
    ax1 = interaction_strength_specific.plot(ax=ax1,
                                             y="frequency",
                                             x="relative_abundance",
                                             kind="scatter",
                                             color=palette[class_name],
                                             label=class_name,
                                             alpha=0.8)
    ax1.set_ylabel("mean frequency",fontsize=12)
    ax1.set_xlabel("Relative abundance (%)",fontsize=12)
    ax1.legend([],[])
    ax1.set_xscale("log")
    ax1.set_ylim([0,1])
    handles1, labels1 = ax1.get_legend_handles_labels()

    ## Process data for violin plot
    combined = pd.concat([interaction_strength_all["frequency"],interaction_strength_specific["frequency"]],axis=1)
    combined = combined.fillna(0)
        
    combined.columns=["all compounds",class_name]
    combined["relative abundance"] = combined.index.map(lambda x: 0 if interaction_strength_all.relative_abundance[x]<10 else 1)
    combined = pd.DataFrame(combined.set_index("relative abundance").stack()).reset_index()
    combined.columns = ["relative abundance","compound class","frequency"]
    
    if get_data:
        return combined
    
    ## Violin plot
    ax2 = sns.violinplot(combined,
                         ax=ax2,
                         x="relative abundance",
                         y="frequency",
                         hue="compound class",
                         cut=0,
                         fill=False,
                         linewidth=1,
                         gap=0.1,
                         palette=palette)
    ax2.set_xticks([0,1],["Below 10%","Above 10%"],fontsize=10)
    ax2.set_ylabel("")
    ax2.set_xlabel("Relative abundance (%)",fontsize=12)
    ax2.legend("")
    ax2.set_ylim([0,1])

    handles2, labels2 = ax2.get_legend_handles_labels()
    handles = [(h1, h2) for h1, h2 in zip(handles1, handles2)]
    return ax1,ax2, handles, labels1