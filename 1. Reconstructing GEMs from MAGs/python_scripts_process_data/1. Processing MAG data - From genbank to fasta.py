#!/usr/bin/env python
# coding: utf-8

# # Change from genebank to fasta-format

# Strategy:
# - Make GEM for all genomes



# ### Make fasta file from genbank

print("**CONVERT FROM GENBANK TO FASTA**")
directory_in_str = "../input/genbank/"
directory = os.fsencode(directory_in_str)
    
gene_ids = {} 

# For each mag


print("\t Converting and saving...")
for file in os.listdir(directory):
    filename = os.fsdecode(file)

    if filename.endswith(".gbk"): 

        gene_id_new = 0
        gene_ids[(filename[:-4],"old_id")] = []
        gene_ids[(filename[:-4],"new_id")] = []
        ofile = open("../output/MAGs_fasta/"+filename[:-4]+".faa","w")
        
        # For each scaffold
        for seq_record in SeqIO.parse("../input/genbank/"+filename,"genbank"):

            # For each gene sequence in the scaffold
            for feature in seq_record.features:

                # If the translation field is in the data -> could be enzyme
                if "translation" in feature.qualifiers.keys():
                    translation = feature.qualifiers["translation"][0]
                    gene_id_old = feature.qualifiers["gene"][0]
                    
                    
                    
                    ofile.write(">"+"gene"+str(gene_id_new)+"\n"+translation+"\n")
                    
                    gene_ids[(filename[:-4],"old_id")].append(gene_id_old)
                    gene_ids[(filename[:-4],"new_id")].append("gene"+str(gene_id_new))
                    
                    gene_id_new = gene_id_new+1
        ofile.close()
gene_ids_df = pd.DataFrame.from_dict({ key:pd.Series(value) for key, value in gene_ids.items()})

gene_ids_df.to_csv("../input/genbank/gene_ids_overview.tsv",sep="\t",index_label=None)

