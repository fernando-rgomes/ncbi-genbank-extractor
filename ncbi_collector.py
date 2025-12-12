from Bio import Entrez, SeqIO
from dotenv import load_dotenv
import pandas as pd
import time
import os

load_dotenv()

Entrez.email = os.getenv("NCBI_EMAIL")

df = pd.read_csv("passeseuarquivo.csv")

#print(df.head())

results = []

for species in df["Especies"]:
    try:
        search = Entrez.esearch(db="nuccore", term=species, retmax=1)
        search_result = Entrez.read(search)

        if not search_result["IdList"]:
            results.append({
                "species": species,
                "locus": None,
                "geo_loc_name": None,
                "collected_by": None,
                "identified_by": None,
                "gene": None
            })
            continue

        seq_id = search_result["IdList"][0]

        handle = Entrez.efetch(db="nuccore", id=seq_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")

        locus = record.name

        # Passo 3 pegar features

        geo = None
        collected = None
        identified = None
        gene = None

        for feat in record.features:
            if feat.type == "source":
                geo = feat.qualifiers.get("geo_loc_name", [None])[0]
                collected = feat.qualifiers.get("collected_by", [None])[0]
                identified = feat.qualifiers.get("identified_by", [None])[0]

            if feat.type == "gene":
                gene = feat.qualifiers.get("gene", [None])[0]

        results.append({
            "species": species,
            "locus": locus,
            "geo_loc_name": geo,
            "collected_by": collected,
            "identified_by": identified,
            "gene": gene
        })

        time.sleep(0.3)




    except Exception as e:
        results.append({
             "species": species,
                "locus": None,
                "geo_loc_name": None,
                "collected_by": None,
                "identified_by": None,
                "gene": None
        })

output = pd.DataFrame(results)
output.to_csv("saida_ncbi.csv")


print("CSV gerado com sucesso!✅")
