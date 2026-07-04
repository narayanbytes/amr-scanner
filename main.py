import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import Bio
from Bio import Seq,SeqIO
import rgi_scanner
import data_cleaner
import os


try:
    organism,file_path = rgi_scanner.loader()
except Exception as e:
    print(f"Setup failed : {e}")
    exit()

# try:
#     data_cleaner.cleaner(flag)
# except Exception as e:
#     print(f"Error : {e}")
#     exit()


try:
    data = pd.read_csv(f"data/{organism}/filtered_data.csv")
except Exception as e:
    print(f"Failed to load data!!\nError: {e}")
    exit()

# print(data.head())["Antibiotic"]["Resistance Mechanism"]["Contig"]

# try:
#     record = SeqIO.read(f"data/{organism}/genome.fna","fasta")
#     print(record.id)
#     print(record.description)
#     genome = record.seq
#     print(f"Length of genome : {len(genome)} base pairs\n")
# except Exception as e:
#     print(f"Failed to load the genome!\nError : {e}")
#     exit()

try:
    genomes = {}
    for record in SeqIO.parse(f"data/{organism}/genome.fna","fasta"):
        genomes[record.id] = record.seq
        print(record.id)
        print(f"Length : {len(genomes[record.id])}\n")
except Exception as e:
    print(f"Failed to load the genome!\nError : {e}")
    exit()   
x_values = data.apply(lambda row: (row["Start"]+row["Stop"])//2,axis = 1).to_list()
y_values = np.ones(len(x_values))
data["DNA Sequence"] = data.apply(lambda row:str(genomes[row["Contig"]][row["Start"]-1 : row["Stop"]-1]),axis = 1)
data["Antibiotic"] = data["Antibiotic"].str.split(";")
data = data.explode("Antibiotic")
data["Antibiotic"] = data["Antibiotic"].str.strip()
Antibiotics = data["Antibiotic"].value_counts().to_frame().reset_index()
Antibiotics.rename(columns = {"count": "Number of resistant genes"},inplace=True)


data["Resistance Mechanism"] = data["Resistance Mechanism"].str.split(";")
data = data.explode("Resistance Mechanism")
data["Resistance Mechanism"] = data["Resistance Mechanism"].str.strip()

# print(data.head()[["Antibiotic","Resistance Mechanism","DNA Sequence"]])
total_ab = data["Antibiotic"].nunique()
print(f"Total antibiotics to which E coli is resistant : {total_ab}")

grouped_data = data.groupby(["Resistance Mechanism"])["Antibiotic"].unique().reset_index()



grouped_data.dropna(subset=["Antibiotic"],inplace=True)
grouped_data["Number of antibiotics"] = grouped_data["Antibiotic"].apply(lambda x: sum((pd.notna(i) and str(i).lower()!="nan") for i in x))
grouped_data["Resistance Mechanism"] = grouped_data["Resistance Mechanism"].str.title()
# print(grouped_data.head())
grouped_data.sort_values(["Number of antibiotics"],ascending = False,inplace=True)
print(grouped_data.head())


# Plotting the data

# contigs = data["Contig"].unique()
# for i,contig in enumerate(contigs):
#     plt.hlines(
#         y = i,xmin = 0,xmax = len(genomes[contig])
#     )
#     temp = data[data["Contig"] == contig]
#     x = (temp["Start"]+temp["Stop"])//2
#     plt.scatter(x,np.full(len(x),i))
#     # plt.yticks(contig)

# plt.yticks(
#     range(len(contigs)),contigs
# )




fig = plt.figure(figsize=(14,10),facecolor = "#ffffff")
gs = GridSpec(2,2)
ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[0,1])
ax3 = fig.add_subplot(gs[1,:])

colors = plt.cm.viridis(np.linspace(0,1,len(grouped_data["Resistance Mechanism"])))
bars = ax1.barh(grouped_data["Resistance Mechanism"],grouped_data["Number of antibiotics"],height = 0.5,color = colors,edgecolor = "black")

for bar in bars:
    # plt.text(x,y,text) --> Put some text at cordinates (x,y)
    ax1.text(bar.get_width(),
    bar.get_y()+bar.get_height()/2,
    str(bar.get_width()),
    color = "Brown",
    fontsize = 16
    )
ax1.set_xlabel("Number Of Antibiotics",fontsize = 12,color = "#212529",fontweight = "bold")
ax1.set_ylabel("Resistance Mechanism",fontsize = 12,color = "#212529",fontweight = "bold")
ax1.set_title("Mechanism Used Against Various Antibiotics",fontweight="heavy",fontsize=16,color = "#2d3748")
ax1.invert_yaxis()
ax1.set_facecolor("#e0f2f1")
# ax1.set_yticks(fontsize = 10)

ax1.grid(axis = "x",linestyle = "--",alpha = 0.4)




# if(len(genomes)==1):
#     ax2.scatter(x_values,y_values,s=50)
#     ax2.set_xlabel("DNA Sequence",fontsize = 12,color = "#212529",fontweight = "bold")
#     ax2.set_title("Position of the resistance genes on the DNA of E Coli",fontweight="bold",color = "#2d3748",fontsize = 16)
#     ax2.set_xlim(0,len(genomes[0].seq))
#     ax2.set_yticks([])
#     ax2.set_facecolor("#e0f2f1")
# else:
contigs = data["Contig"].unique()
for i,contig in enumerate(contigs):
    ax2.hlines(
        y = i,xmin = 0,xmax = len(genomes[contig])
    )
    temp = data[data["Contig"] == contig]
    x = (temp["Start"]+temp["Stop"])//2
    ax2.scatter(x,np.full(len(x),i))
    # plt.yticks(contig)

ax2.set_yticks(
    range(len(contigs)),contigs
)
ax2.set_xlabel("Length of the various contigs",fontsize=12,color = "#212529",fontweight = "bold")
ax2.set_title(f"Position of the resistance genes on the genome of\n{organism.replace("_"," ")}",fontweight="bold",color = "#2d3748",fontsize = 16)
ax2.set_ylabel("Contigs",fontsize = 12,color = "#212529",fontweight = "bold")

bars = ax3.bar(Antibiotics["Antibiotic"],Antibiotics["Number of resistant genes"],edgecolor = "black")


for bar, antibiotic in zip(bars,Antibiotics["Antibiotic"]):
    ax3.text(bar.get_x() + bar.get_width()/2,
    bar.get_height()/2,
    antibiotic,
    ha = "center",
    va = "center",
    rotation = 90,
    color = "black",
    fontweight="normal")

for bar in bars:
    ax3.text(bar.get_x()+bar.get_width()/2,
    bar.get_height(),
    str(bar.get_height()),
    color = "brown",
    ha = "center"
    )
ax3.set_xlabel("Antibiotics",fontsize = 12,color = "#212529",fontweight = "bold")
ax3.set_ylabel("Number of resistant genes\n(conferring resistance to a specific antibiotic)",fontsize = 12,color = "#212529",fontweight = "bold")
ax3.set_title("Resistance Genes in E coli against various antibiotics",fontweight = "bold",fontsize = 16, color = "#2d3748")
ax3.set_xticks([])
ax3.grid(axis = "y",linestyle = "--",alpha = 0.4)
ax3.set_facecolor("#e0f2f1")
ax3.text(
    0.98,0.98,
    f"{organism.replace("_"," ")}\nTotal Antibiotics to which bacteria is resitant : {total_ab}",
    transform = ax3.transAxes,
    ha = "center",
    va = "top",
    bbox = dict(
        facecolor = "white",
        edgecolor = "black",
        alpha = 0.8
    )
)
fig.tight_layout(pad = 2)
plt.subplots_adjust(left = 0.2)
plt.show()