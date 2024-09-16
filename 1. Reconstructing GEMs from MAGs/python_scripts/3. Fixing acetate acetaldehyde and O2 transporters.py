print("**FIXING TRANSPORT REACTIONS FOR ACETALDEHYDE, ACETATE, AND OXYGEN**")

with open("../output/relevant_MAGs_99.txt") as text_file:
    relevant_MAGs = text_file.read().split("\n")

relevant_MAGs = [string.replace("\t","") for string in relevant_MAGs]


##### Load models
print("\t Load models...")

GEMs_dict = {}

directory = os.fsencode("../output/GEMs/GEMs_intermediate/GEMs_adapt/")

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    
    if filename.endswith(".xml") and filename[:-4] in relevant_MAGs:
        GEMs_dict[filename[:-4]]= reframed.load_cbmodel("../output/GEMs/GEMs_intermediate/GEMs_adapt/"+filename)


##### Find acetate producers/consumers

print("\t Process models...")

#  Find which models has GPRs for producing acetate through phosphotransacetylase(R_PTAr) and acetate kinase (R_ACKr)

has_enzymes_for_acetate = {}
for MAG in GEMs_dict.keys():
    
    has_enzyme = []
    if "R_ACKr" in GEMs_dict[MAG].reactions:
        if GEMs_dict[MAG].reactions["R_ACKr"].gpr!=None:
            has_enzyme.append("R_ACKr")
    if "R_PTAr" in GEMs_dict[MAG].reactions:
        
        if GEMs_dict[MAG].reactions["R_PTAr"].gpr!=None:
            has_enzyme.append("R_PTAr")
    
    has_enzymes_for_acetate[MAG]= len(has_enzyme)==2

has_enzymes_for_acetate_MAGs = list(pd.Series(has_enzymes_for_acetate)[pd.Series(has_enzymes_for_acetate)].index)

##### Find hits on the ACt2r protein from the TCDB database
print("\t Process blast results...")

ACt2r_MAGs = []

ACt2r_MAGs_data = []
for filename in os.listdir("transporters/"):
    
    if filename.endswith(".tsv"):

        transport= pd.read_csv("../output/transporters/"+filename,sep="\t",header=None)
        transport.columns = ["query acc.ver", "subject acc.ver", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score"]
        transport.reset_index(inplace=True)
        
        if transport.shape[0]>1:
            
            transport.sort_values(by="evalue",ascending=True,inplace=True)
            transport.reset_index(inplace=True)
            
            query = transport.loc[0,"query acc.ver"]
            gene = transport.loc[0,"subject acc.ver"]
            best_evalue = transport.loc[0,"evalue"]
            bit_score = transport.loc[0,"bit score"]
            
            if best_evalue<1e-5 and bit_score>20:
                if "2.A.1.13.1" in query:
                    ACt2r_MAGs.append(filename[:-4])
                    
                    ACt2r_MAGs_data.append((filename[:-4],gene,best_evalue,bit_score,"2.A.1.13.1"))

                elif "2.A.21.7.3" in query: 
                    ACt2r_MAGs.append(filename[:-4])
                    ACt2r_MAGs_data.append((filename[:-4],gene,best_evalue,bit_score,"2.A.21.7.3"))


ACt2r_MAGs_df = pd.DataFrame(ACt2r_MAGs_data,columns=["MAG","gene","evalue","bit_score","TCDB_id"]).sort_values("evalue")
ACt2r_MAGs_df.set_index("MAG",inplace=True)

##### Add acetate transport

# Two different conditions.

# 1. The MAG has a hit in the TCDB database -> add reversible reaction
# 2. The MAG has a hit for the enzymes in acetate production -> add only producing reaction

print("\t Change models...")

GEMs_dict2 = {}

for MAG,model in GEMs_dict.items():
    model_copy = model.copy()
    
    if "R_Acabc" in model_copy.reactions.keys(): 
        model_copy.remove_reaction("R_Acabc")
      
    if "R_ACt2r" in model_copy.reactions.keys():
        model_copy.remove_reaction("R_ACt2r")
        
    if MAG in ACt2r_MAGs_df.index.values:
        model_copy.add_reaction_from_str("R_ACt2r: M_ac_c + M_h_c --> M_ac_e + M_h_e")
        model_copy.reversible=True
        model_copy.lb=-1000
        
        GPR = parse_gpr_rule(ACt2r_MAGs_df.loc[MAG,"gene"])
        model_copy.set_gpr_association("R_ACt2r",GPR)
        
        
    elif has_enzymes_for_acetate_MAGs:
        model_copy.add_reaction_from_str("R_ACt2r: M_ac_c + M_h_c --> M_ac_e + M_h_e")
        
    GEMs_dict2[MAG]=model_copy

### Removing oxygen related reactions

GEMs_dict3 = {}
for MAG, model in GEMs_dict2.items():
    if "M_o2_e" not in model.metabolites:
        continue
        
    model_copy = model.copy()
    model_copy.remove_reactions(model_copy.get_metabolite_reactions("M_o2_e"))
    model_copy.remove_metabolite("M_o2_e")
    model_copy.remove_reaction("R_EX_o2_e")
    
    GEMs_dict3[MAG]=model_copy

##### Save models
print("\t Save models...")

for MAG,model in GEMs_dict3.items():
    model.update()
    reframed.save_cbmodel(model,"../output/GEMs/GEMs_intermediate/GEMs_ACt2r/"+MAG+".xml")