from Bio import Entrez, SeqIO
from dotenv import load_dotenv
import pandas as pd
import time
import os

load_dotenv()

Entrez.email = os.getenv("NCBI_EMAIL")

df = pd.read_csv("passeaquioarquivo.csv")

print(df.head())

results = []

'''for species in df["Especies"]:
    try:
        search = Entrez.esearch(db="nuccore", term=species, retmax=1)
        search_result = Entrez.read(search)

        if not search_results["IdList"]:
            resultados.append({
                "species": especie,
                "locus": None,
                "geo_loc_name": None,
                "collected_by": None,
                "identified_by": None,
                "gene": None
            })
            continue

        seq_id = search_result["IdList"][0]


    except'''


    
