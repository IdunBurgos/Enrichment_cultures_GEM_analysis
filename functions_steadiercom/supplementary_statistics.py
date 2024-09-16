
import pandas as pd
import scipy.stats as stats
import math


def find_non_dependent(row,metric,family_groups):
    return len(family_groups[row.name[0]]) - row[metric] #row.flux_mg


def statistics_adjustments(statistics_df):
    
    statistics_df = statistics_df.sort_values(by="p_value").copy()
    statistics_df["i"] = statistics_df["p_value"].rank(method="max")
    statistics_df["p_value_benjamini_h"] = statistics_df.apply(lambda row: min(row.p_value*statistics_df.shape[0]/row.i,1),axis=1)
    statistics_df.sort_index(inplace=True)
    return statistics_df


def statistics_function(steadier_sample_cross,dependent_variable,independent_variable,family_groups,metric="flux_mg",metric_thresh=1e-6,pvalue_thresh=0.1):

    # Get average of each family according to each possible value of the independent variable
    # dependent_variable,dependent_variable.split("_")[1] here it decides if it is in the receiver or in the donor (dependent_variable.split("_")[1]) and groups by the [family_receiver,receiver,compound] and takes the mean of this
    steadiercom_crossfeeding_donor = steadier_sample_cross.loc[:,[dependent_variable,dependent_variable.split("_")[1],independent_variable,metric]].groupby([dependent_variable,dependent_variable.split("_")[1],independent_variable]).mean().copy()
    dependent = steadiercom_crossfeeding_donor[steadiercom_crossfeeding_donor[metric]>metric_thresh].reset_index().groupby([dependent_variable,independent_variable]).count().copy()
    not_dependent = dependent.apply(find_non_dependent,metric=metric,family_groups=family_groups,axis=1)

    # Add data for the missing values
    all_categories =set(not_dependent.index.get_level_values(1))

    for family in dependent.index.get_level_values(0):
        for category in all_categories-set(not_dependent.xs(family).index):
            not_dependent[(family,category)]=len(family_groups[family]) 
            
    concat_df = pd.concat({"dependent":dependent[metric],"not_dependent":not_dependent},axis=1).fillna(0)


    statistics = {}
    for independent_var in set(concat_df.index.get_level_values(1)):
        # Get the sub_df for the super class
        concat_df_sub = concat_df.xs(independent_var,level=1).copy()

        statistics[independent_var] = {}

        # For each row (each family)
        for i,row in concat_df_sub.iterrows():
            
            # Get the data for all other family
            other = concat_df_sub.loc[concat_df_sub.index[concat_df_sub.index!=i],:]
            
            data = pd.DataFrame({i:row,"other":other.sum()}).transpose().to_numpy()
            statistics[independent_var][(i,"data")]= data
            
            
            # Get odds ratio
            odds_ratio_num = data[0][0]
            odds_ratio_den = data[0][0] + data[0][1]
            other_num = data[1][0]
            other_den = data[1][0] + data[1][1]
            
            if odds_ratio_den == 0 or other_den == 0 or other_num==0:
                odds_ratio = math.inf
            else:
                odds_ratio = (odds_ratio_num / odds_ratio_den) / (other_num / other_den)
            statistics[independent_var][(i, "odds_ratio")] = odds_ratio
            
            # Calculate the Barnard exact statistical 
            p_value = stats.barnard_exact(pd.DataFrame({i:row,"other":other.sum()}).transpose().to_numpy(),alternative="greater")
            statistics[independent_var][(i,"p_value")]= p_value.pvalue
            

    statistics_df = pd.DataFrame(statistics)
    category_values = statistics_df.xs('p_value', level=1)
    

    values = {}
    for family,independent_var in category_values[category_values[category_values<pvalue_thresh].notnull()].stack().index:
        values[family,independent_var]= {"p_value":statistics_df.loc[(family,"p_value"),independent_var],"odds_ratio":statistics_df.loc[(family,"odds_ratio"),independent_var],"data":statistics_df.loc[(family,"data"),independent_var],"# interaction":concat_df.loc[(family,independent_var),"dependent"],"# no interaction":concat_df.loc[(family,independent_var),"not_dependent"]}

    significant = pd.DataFrame(values).transpose()
    if significant.empty==False:
        significant.index.names = (dependent_variable,independent_variable)

    return significant