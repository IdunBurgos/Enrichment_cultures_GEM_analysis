import reframed

import pandas as pd

from sklearn.metrics.pairwise import pairwise_distances
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc
import seaborn as sns

from matplotlib.pyplot import gcf
import matplotlib.pyplot as plt 

import functions.colors_MAGs
import functions.general_functions as general_func


chebi_lut,chebi_interesting,chebi_colors_ser = functions.colors_MAGs.chebi_rxn_color_func()

def MAG_MAG_clustermap_data(GEMs_dict,only_exchange_rxns=False):
    """
    Args:
        GEMs_dict (dict): A dictionary where keys are model IDs and values are GEMs.
        only_exchange_rxns (bool): If True, only considers exchange reactions; otherwise, considers all reactions.
        
        
    Returns:
        jaccard (pd.DataFrame): A DataFrame representing the Jaccard distance matrix.
        linkage (ndarray): An array representing the hierarchical clustering linkage.
        rxns_models (list): A list of lists, each containing the reactions of a model.
    
    """
    
    if only_exchange_rxns:
        rxns_models = [list(model.get_exchange_reactions()) for model in GEMs_dict.values()]
    else:
        rxns_models = [list(model.reactions.keys()) for model in GEMs_dict.values()]
    
    
    assert len(rxns_models)==len(GEMs_dict.keys()),"You are not processing data for all models. Check that loading was successful."
        
    rxns_union = set().union(*rxns_models)
    rxns_dict = {model_id:{rxn:rxn in model.reactions for rxn in rxns_union} for model_id, model in GEMs_dict.items()}

    df = pd.DataFrame(rxns_dict)
    
    assert df.shape == (len(rxns_union), len(GEMs_dict)), "DataFrame dimensions are incorrect."
    
    # Compute distance and get dataframe
    jac_dist = pairwise_distances(df.T.to_numpy(), metric='jaccard')

    # Make this as dataframe
    jaccard = pd.DataFrame(jac_dist,index=df.columns,columns=df.columns)
    
    assert jac_dist.shape == (len(GEMs_dict), len(GEMs_dict)), "Jaccard distance matrix dimensions are incorrect."
    
    # Calculate linkage
    linkage = hc.linkage(sp.distance.squareform(jaccard), method='average')
    
    return jaccard,linkage



def MAG_rxn_clustermap_data(GEMs_dict,interesting_super_classes=False):
    
    if interesting_super_classes:
        chebi_interesting = general_func.chebi_selected(interesting_super_classes=interesting_super_classes)
    else:
        chebi_interesting = general_func.chebi_selected()
        
    relevant_rxns = chebi_interesting.index
    
    rxns_models = [list(model.get_exchange_reactions()) for model in GEMs_dict.values()]

    rxns_union = set().union(*rxns_models)

    rxns_dict = {model_id:{rxn:rxn in model.reactions for rxn in rxns_union} for model_id, model in GEMs_dict.items()}

    df = pd.DataFrame(rxns_dict)
    
    
    df=df.loc[df.sum(axis=1)<df.shape[1]*0.95,:].copy()
    df = df.loc[df.sum(axis=1)>df.shape[1]*0.05,:].copy()
    df = df.loc[df.index.str.contains("EX_"),:].copy()

    # only get reactions of relevant classes
    df = df.loc[df.index.isin(relevant_rxns),:].copy()


    # Distance: This finds the distance between the models
    jac_dist_mags = pairwise_distances(df.T.to_numpy(), metric='jaccard')
    jaccard_mags = pd.DataFrame(jac_dist_mags,index=df.columns,columns=df.columns)

    jac_dist_rxns = pairwise_distances(df.to_numpy(), metric='jaccard')
    jaccard_rxns = pd.DataFrame(jac_dist_rxns,index=df.index,columns=df.index)

    # Linkage
    linkage_row = hc.linkage(sp.distance.squareform(jaccard_rxns), method='average')
    linkage_col = hc.linkage(sp.distance.squareform(jaccard_mags), method='average')
    
    return df,jaccard_mags,jaccard_rxns,linkage_row,linkage_col


def producers_consumers_sim(GEMs_dict):

    GEMS_consumers = {}
    GEMs_producers = {}

    for GEM_id,model in GEMs_dict.items():
        complete_env = reframed.Environment.from_model(model)

        sol = reframed.FVA(model,constraints=complete_env,reactions=model.get_exchange_reactions())

        GEMS_consumers[GEM_id] = {}
        GEMs_producers[GEM_id] = {}
        print(GEM_id)
        for rxn,sol_rxn in sol.items():
            if sol_rxn[0]<-1e-6:
                GEMS_consumers[GEM_id][rxn]=sol_rxn[0]
            else:
                GEMS_consumers[GEM_id][rxn]=0

            if sol_rxn[1]>1e-6:
                GEMs_producers[GEM_id][rxn]=sol_rxn[1]
            else:
                GEMs_producers[GEM_id][rxn]=0
    
    with pd.option_context("future.no_silent_downcasting", True):
        GEMs_producers_df = pd.DataFrame(GEMs_producers).fillna(0).infer_objects(copy=False)
    
    with pd.option_context("future.no_silent_downcasting", True):
        GEMS_consumers_df = pd.DataFrame(GEMS_consumers).fillna(0).infer_objects(copy=False)
    
    
    
    return GEMS_consumers_df,GEMs_producers_df

"""
def MAGs_rxns_data_processing_sim(GEMS_df,relevant_rxns=chebi_colors_ser.index):

    
    
    GEMS_df = GEMS_df.apply(abs)
    # Filtering
    
    
    GEMS_df_copy = GEMS_df.loc[GEMS_df.index.isin(relevant_rxns),:].copy()
    
    
    #df,jaccard_mags,jaccard_rxns,linkage_row,linkage_col = MAG_rxn_clustermap_data(GEMs_dict)
    
    #GEMS_df_copy.loc[GEMS_df_copy.index.isin(df.index)].copy()
    
    
    GEMS_df_copy = GEMS_df_copy.loc[GEMS_df_copy.sum(axis=1)>(GEMS_df_copy.shape[1]*0.05*10),:].copy()
    #GEMS_df_copy = GEMS_df_copy.loc[GEMS_df_copy.sum(axis=1)<(GEMS_df_copy.shape[1]*0.95*10),:].copy()
    
    # Distance 
    dist_mag = pairwise_distances(GEMS_df_copy.T.to_numpy(), metric='euclidean')
    distance_mags = pd.DataFrame(dist_mag,index=GEMS_df_copy.columns,columns=GEMS_df_copy.columns)
    assert distance_mags.shape[0]==distance_mags.shape[1], "jaccard_mags: "+str(distance_mags.shape) 
    
    distance_mags = (distance_mags.T+distance_mags)/2
    
    assert distance_mags.shape[0]==GEMS_df.shape[1] & distance_mags.shape[1]==GEMS_df.shape[1], "distance_mags: " + str(distance_mags.shape)+", GEMS_df: " + str(GEMS_df.shape)
    
    dist_rxns = pairwise_distances(GEMS_df_copy.to_numpy(), metric='euclidean')
    distance_rxns = pd.DataFrame(dist_rxns,index=GEMS_df_copy.index,columns=GEMS_df_copy.index)
    assert distance_rxns.shape[0]==distance_rxns.shape[1]
    
    distance_rxns = (distance_rxns.T+distance_rxns)/2
    assert distance_rxns.shape[0]==GEMS_df_copy.shape[0] & distance_rxns.shape[1]==GEMS_df_copy.shape[0],"distance_rxns: " + str(distance_rxns.shape)+", GEMS_df_copy: " + str(GEMS_df_copy.shape)
    
    # Linkage
    linkage_row = hc.linkage(sp.distance.squareform(distance_rxns), method='ward')
    linkage_col = hc.linkage(sp.distance.squareform(distance_mags), method='ward')
    
    return GEMS_df_copy, distance_mags,distance_rxns,linkage_col,linkage_row



"""

def MAGs_rxns_data_processing_sim(GEMS_df,interesting_super_classes=False):
    """
    Args:
        GEMS_df (pd.DataFrame): Dataframe with values of either the consumed or produced compounds
        
        
    Returns:
        GEMS_df_copy
        distance_mags (pd.DataFrame)
        distance_rxns (pd.DataFrame)
        linkage_col (ndarray)
        linkage_row (ndarray)
    """
    if interesting_super_classes:
        chebi_interesting = general_func.chebi_selected(interesting_super_classes=interesting_super_classes)
    else:
        chebi_interesting = general_func.chebi_selected()
        
    relevant_rxns = chebi_interesting.index
    
    GEMS_df = GEMS_df.apply(abs)
    
    # Filtering
    GEMS_df_copy = GEMS_df.loc[GEMS_df.index.isin(relevant_rxns),:].copy()
    
    GEMS_df_copy = GEMS_df_copy.loc[GEMS_df_copy.sum(axis=1)>(GEMS_df_copy.shape[1]*0.05*10),:].copy()
     
    return GEMS_df_copy






    
    