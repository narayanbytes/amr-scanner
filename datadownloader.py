from Bio import Entrez
from Bio import Seq, SeqIO
import pandas as pd
import requests
import gzip
import shutil
import os
Entrez.tool = "Antibiotic Reistance Analyzer"
def downloader():
    name = input("Enter the name of the bacteria : ").replace(" ","_").title()
    file_path = f"data/{name}"
    if os.path.exists(f"{file_path}"):
        return name,None
        
    try:
        print("\nThe data for the organism needs to be downloaded!\n")
        email = input("Please enter your email :  ").lower()
        try:
            os.makedirs(f"{file_path}",exist_ok = True)
            Entrez.email = email
            handle = Entrez.esearch(
                db = "assembly",
                term = f"{name}[Organism]",
                retmax = 5
            )

        except Exception as e:
            print(f"Error : {e}")
            exit()

        data = Entrez.read(handle)
        handle = Entrez.esummary(db = "assembly",id = data["IdList"])
        final_data = Entrez.read(handle)
        # print(final_data["DocumentSummarySet"]["DocumentSummary"][0].keys())
        print("Please select an organism (with a valid FTP Path) : ")
        for i,assembly in enumerate(final_data["DocumentSummarySet"]["DocumentSummary"],start = 1):
            print(f"{i}-->\n Name: {assembly["Organism"]}\n Assembly: {assembly["AssemblyAccession"]}\n FTP Path : {assembly["FtpPath_RefSeq"]}\n FTP Path2 : {assembly["FtpPath_GenBank"]}")
        print()
        res = int(input("Enter your response here : (Select an option out of 1-5):  "))
        while res not in range(1,6):
            print("Invalid Response!")
            res = int(input("Enter a valid response : "))
        ftp_path1 = final_data["DocumentSummarySet"]["DocumentSummary"][res-1]["FtpPath_RefSeq"]
        ftp_path2 = final_data["DocumentSummarySet"]["DocumentSummary"][res-1]["FtpPath_GenBank"]
        

        if ftp_path1 != "":
            ftp_path = ftp_path1
        elif ftp_path2 != "":
            ftp_path = ftp_path2
        else:
            print("No downloadable data exists! Please try another organism")
            exit()
        

        ftp_path = ftp_path.replace("ftp://","https://")
        folder = ftp_path.split("/")[-1]
        url = f"{ftp_path}/{folder}_genomic.fna.gz"
        response = requests.get(url)

        with open(f"{file_path}/genome.fna.gz","wb") as f:
            f.write(response.content)

        compressed_file = gzip.open(f"{file_path}/genome.fna.gz","rb")

        normal = open(f"{file_path}/genome.fna","wb")
        shutil.copyfileobj(compressed_file,normal)
        compressed_file.close()
        normal.close()
        os.remove(f"{file_path}/genome.fna.gz")
        # organism = final_data["DocumentSummarySet"]["DocumentSummary"][res-1]["Organism"]
        organism = name
        return organism,file_path
        # for record in SeqIO.parse(f"{file_path}/genome.fna","fasta"):
        #     print(record.id)
        #     print(record.description)
        #     genome = record.seq
        #     print(len(genome))
    except Exception as e:
        os.rmdir(file_path)
        print(f"An error occcured:{e}")
        exit()
