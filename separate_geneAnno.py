import pandas as pd
file_in = "PATH/PGCNA_work/PostPGCNA_normB/MEV_Vis/geneModuleWeightedDegree_normB_0.3_GEO_ModCon.txt" # change path
df = pd.read_csv(file_in, sep="\t")

print(df.head())
print(df.columns)
print(df["Gene"].head())

left_side = [v.split("#")[0] for v in df["Gene"].values]
right_side = [v.split("#")[1] for v in df["Gene"].values]

df["Gene_1"] = left_side
df["ID"] = right_side
file_out = file_in.replace(".txt", "_out.txt")
df.to_csv(file_out, sep="\t",index=False)
