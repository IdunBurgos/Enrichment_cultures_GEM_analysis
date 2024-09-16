import reframed
import pandas as pd
import copy

def enrich_syncon_environments():
    """
    Creates syncon environments: SC1_X, SC2_X, SC1_C, and SC2_X
    
    ---
    
    From SynCon2 to SynCon1

    SynCon2 had the same composition, but it lacked ascorbic acid, and it had two additional trace elements (selenium and wolfram)

    Add
    - Ascorbate: ascb__L

    Remove
    - Selenite (inorganic selenium): slnt
    - (wolfram: not in BiGG database)
    
    """
    
    substrate_composition = {"Cellulose":["cellb","cell3","cell4","cell5"],
                             "Xylan":["xylb","xyl3","xylan4","xylan8"]} 


    syncon2 = list(pd.read_csv("../input/syncon2media_combined.csv",header=None)[0].values)

    syncon1 = copy.copy(syncon2)
    syncon1.remove("slnt")
    syncon1.append("ascb__L")

    media = {}

    for polysaccharide,compounds in substrate_composition.items():
        media["SC1_"+polysaccharide[0]] = copy.copy(syncon1)
        media["SC1_"+polysaccharide[0]].extend(compounds)

        media["SC2_"+polysaccharide[0]] = copy.copy(syncon2)
        media["SC2_"+polysaccharide[0]].extend(compounds)
        
    return media


def community_syncon_environments():
    
    """
    Assigns the communities to the syncon environment used in their media.
    
    """
    media = enrich_syncon_environments()
    
    MAG2community_id = pd.read_csv("../output/MAG2community_id.tsv",sep="\t",header=None)
    communities = MAG2community_id[1].unique()
    syncon_environments = {}

    for community_id in communities:
        if community_id in  ["CM_A","CM_P","M_A","M_P"]:
            syncon_environments[community_id]=reframed.Environment.from_compounds(media["SC2_C"])
        elif community_id in ["CD_A","CD_P"]:
            syncon_environments[community_id]=reframed.Environment.from_compounds(media["SC1_C"])
        elif community_id.endswith("_X"):
            syncon_environments[community_id] = reframed.Environment.from_compounds(media["SC1_X"])
    return syncon_environments

