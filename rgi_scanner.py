import pandas as pd
import subprocess
from datadownloader import downloader
from data_cleaner import cleaner
def loader():
    try:
        organism,file_path = downloader()
        
        if file_path == None:
            return organism,None
    except Exception as e:
        print(f"Error while downloading the data! : {e}")
        exit()
    input_genome = f"{file_path}/genome.fna"
    output_data = f"{file_path}/pipeline_results"

    # Original rgi CLI comamand --> rgi main --input_sequence input_genome --output_file output_data --local --clean

    rgi_command = [
        "rgi","main",
        "--input_sequence",input_genome,
        "--output_file",output_data,
        # "--local",
        "--clean",
        "--input_type","contig"
    ]

    print(f"<<<--------Starting the RGI engine for {organism.replace("_"," ")}-------->>>\n")
    print("Please wait....It may take a few minutes.")
    process = subprocess.run(rgi_command,capture_output=True,text=True)

    if(process.returncode == 0):
        print("Execution Successful!")
        cleaner(file_path)
        return organism,file_path
        # print(process.stdout)
        # print(process.stderr)
        # data = pd.read_csv(f"{output_data}.txt",sep="\t")

        # print(data.head(10))
        # print(data.columns)
    else:
        print("Error!")
        print(process.stderr)