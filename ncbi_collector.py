from Bio import Entrez, SeqIO
from dotenv import load_dotenv
from tqdm import tqdm
import pandas as pd
import time
import os

# ===============================
# CONFIGURAÇÃO
# ===============================

load_dotenv()
Entrez.email = os.getenv("NCBI_EMAIL")

df = pd.read_csv("species/chile.csv")

genes_mitocondriais = [
    "nad1", "nad2", "nad3", "nad4", "nad4L", "nad5", "nad6",
    "cob", "cytb",
    "cox1", "COI", "cox2", "cox3",
    "atp6", "atp8",
    "16S", "16S rRNA", "16s",
    "12S", "12S rRNA", "12s",
    "trnA", "trnC", "trnD", "trnE", "trnF", "trnG", "trnH",
    "trnI", "trnK", "trnM", "trnN", "trnP", "trnQ", "trnR",
    "trnT", "trnV", "trnW", "trnY",
    "trnL1", "trnL2",
    "trnS1", "trnS2"
]

results = []

# ===============================
# FUNÇÃO DE SCORE
# ===============================

def score_registro(geo, collected, identified):
    return sum(v is not None for v in [geo, collected, identified])

# ===============================
# LOOP PRINCIPAL
# ===============================

for species in tqdm(df["Especies"], desc="Espécies", position=0):

    # 🔎 Buscar TaxID
    tax_search = Entrez.esearch(
        db="taxonomy",
        term=f"{species}[Scientific Name]"
    )
    tax_result = Entrez.read(tax_search)

    taxid = tax_result["IdList"][0] if tax_result["IdList"] else None

    encontrou_algo = False

    # 🔁 Barra interna para genes
    for gene in tqdm(
        genes_mitocondriais,
        desc=f"Genes ({species})",
        leave=False,
        position=1
    ):

        # 🔁 Montar termo de busca
        if taxid:
            termo_busca = f"txid{taxid} AND {gene}"
        else:
            termo_busca = f"{species}[Organism] AND {gene}"

        try:
            search = Entrez.esearch(
                db="nuccore",
                term=termo_busca,
                retmax=200
            )
            search_result = Entrez.read(search)

            if not search_result["IdList"]:
                continue

            encontrou_algo = True
            melhor_score = -1
            melhor_registro = None

            # 🔎 Avaliar apenas os primeiros IDs
            for seq_id in search_result["IdList"][:20]: # Limitando pesquisa

                handle = Entrez.efetch(
                    db="nuccore",
                    id=seq_id,
                    rettype="gb",
                    retmode="text"
                )
                record = SeqIO.read(handle, "genbank")

                locus = record.name
                geo = collected = identified = None

                for feat in record.features:
                    if feat.type == "source":
                        geo = feat.qualifiers.get("geo_loc_name", [None])[0]
                        collected = feat.qualifiers.get("collected_by", [None])[0]
                        identified = feat.qualifiers.get("identified_by", [None])[0]

                score = score_registro(geo, collected, identified)

                if score > melhor_score:
                    melhor_score = score
                    melhor_registro = {
                        "species": species,
                        "taxid": taxid,
                        "gene": gene,
                        "busca": termo_busca,
                        "n_sequencias": len(search_result["IdList"]),
                        "locus_representativo": locus,
                        "geo_loc_name": geo,
                        "collected_by": collected,
                        "identified_by": identified
                    }

                # ⚡ Atalho: registro perfeito
                if score == 3:
                    break

                time.sleep(0.2)

            if melhor_registro:
                results.append(melhor_registro)

            time.sleep(0.4)

        except Exception:
            continue

    # ❗ Espécie sem nenhum resultado
    if not encontrou_algo:
        results.append({
            "species": species,
            "taxid": taxid,
            "gene": None,
            "busca": None,
            "n_sequencias": 0,
            "locus_representativo": None,
            "geo_loc_name": None,
            "collected_by": None,
            "identified_by": None
        })

# ===============================
# SAÍDA FINAL
# ===============================

output = pd.DataFrame(results)
output.to_csv("saida_ncbi_mitocondrial.csv", index=False)

print("CSV mitocondrial gerado com sucesso! ✅")
