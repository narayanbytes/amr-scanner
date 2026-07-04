import pandas as pd
import numpy as numpy

def cleaner(file_path):
    try:
        data = pd.read_csv(f"{file_path}/pipeline_results.txt",sep="\t")
        data = data.drop(columns=["AST_Source","Orientation","SNPs_in_Best_Hit_ARO","Other_SNPs","Predicted_DNA", "Predicted_Protein","CARD_Protein_Sequence","Hit_Start","Hit_End","Nudged","Note","ID","Model_ID"])
        data.drop(data[data["Cut_Off"]=="Loose"].index,inplace = True)
        data.rename(columns={"Percentage Length of Reference Sequence":"% Length of Reference Sequence"},inplace=True)
        # print(data["Start"].dtype)
        # print(len(data))
        # print(data.head()[['ORF_ID', 'Contig', 'Cut_Off',  'Best_Hit_Bitscore', 'Best_Hit_ARO', 'Best_Identities', 'Model_type', 'Drug Class', 'AMR Gene Family', '% Length of Reference Sequence', 'Antibiotic']])
        data.to_csv(f"{file_path}/filtered_data.csv",index=False)
    except Exception as e:
        print(f"An error occured while cleaning the file:\n{e}")