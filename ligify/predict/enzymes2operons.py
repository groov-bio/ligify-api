import json
import os
import re
import time
import requests
from functools import lru_cache

from utils import make_request  # Assuming make_request uses a persistent session and rate limiting

# Create a global session for all external requests
session = requests.Session()
session.headers.update({"User-Agent": "my-bioinfo-app/1.0"})

# If your make_request wrapper now returns a requests.Response but uses the session internally,
# you might want to unify all external calls through it for uniform timing and reuse.
# For example, replace direct requests.get() calls with make_request if it uses `session` internally.

@lru_cache(maxsize=None)
def fetch_uniprot_json(query_url: str):
    # A cached fetch function for Uniprot responses
    resp = session.get(query_url)
    resp.raise_for_status()
    return resp.json()

@lru_cache(maxsize=None)
def fetch_protein_fasta(accession: str):
    # Use make_request for NCBI calls since you have rate limit logic there
    URL = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="
        + accession
        + "&rettype=fasta"
        + f"&api_key={os.getenv('NcbiApiKey')}"
    )
    resp = make_request("GET", URL)  # ensure make_request uses the shared session if possible
    resp.raise_for_status()
    seq = "".join(line for line in resp.text.split("\n")[1:])
    return seq

def protein2chemicals(accession: str):
    url = "https://rest.uniprot.org/uniprotkb/search?query=" + accession + "&format=json"
    try:
        protein = fetch_uniprot_json(url)
    except Exception:
        return None

    if len(protein.get("results", [])) > 0:
        protein_comments = protein["results"][0].get("comments", [])
        protein_data = {}

        # Extract various data from the comments
        response_function = [c["texts"][0]["value"] for c in protein_comments if c["commentType"] == "FUNCTION"]
        if response_function:
            protein_data["function"] = response_function[0]

        catalysis = [c["reaction"]["name"] for c in protein_comments if c["commentType"] == "CATALYTIC ACTIVITY"]
        if catalysis:
            protein_data["catalysis"] = catalysis[0]

        # Attempt to fetch ligands
        try:
            reaction = [
                c["reaction"]["reactionCrossReferences"] 
                for c in protein_comments if c["commentType"] == "CATALYTIC ACTIVITY"
            ]
            if reaction:
                LIGANDS = [r["id"] for r in reaction[0] if r["database"] == "ChEBI"]
                if LIGANDS:
                    protein_data["ligands"] = LIGANDS
        except Exception:
            pass

        induction = [c["texts"][0]["value"] for c in protein_comments if c["commentType"] == "INDUCTION"]
        if induction:
            protein_data["induction"] = induction[0]

        pathway = [c["texts"][0]["value"] for c in protein_comments if c["commentType"] == "PATHWAY"]
        if pathway:
            protein_data["pathway"] = pathway[0]

        return protein_data
    return None

@lru_cache(maxsize=None)
def fetch_uniprot_reg_data(accession: str):
    url = "https://rest.uniprot.org/uniprotkb/search?query=" + accession + "&format=json"
    try:
        data = fetch_uniprot_json(url)["results"][0]
        dois = []
        for j in data.get("references", []):
            doi = None
            title = None
            if "citationCrossReferences" in j["citation"]:
                for k in j["citation"]["citationCrossReferences"]:
                    if k["database"] == "DOI":
                        doi = k["id"]
                        title = j["citation"]["title"]
                        break
                if doi and title:
                    dois.append({"doi": doi, "title": title})

        regulator = {
            "annotation": data["features"][0]["description"] if data.get("features") else "No data available",
            "id": data["primaryAccession"],
            "references": dois if dois else "No data available",
            "length": data["sequence"]["length"] if data.get("sequence") else "No data available",
        }
    except Exception:
        regulator = {
            "annotation": "No data available",
            "id": "No data available",
            "references": "No data available",
            "length": "No data available",
        }

    return regulator

def pull_regulators(protein, rxn):
    regulator_pattern = re.compile(r"regulator|repressor|activator")

    reg_data = []
    ligand_names = []

    # Cache or store operon calls if same operon is processed multiple times
    if protein.get("context", "EMPTY") != "EMPTY":
        operon = protein["context"]["operon"]
        # Pre-fetch all sequences and data for operon genes if possible
        # This reduces repeated calls in a loop

        for gene in operon:
            if "description" in gene and regulator_pattern.search(gene["description"]):
                entry = {
                    "refseq": gene["accession"],
                    "annotation": gene["description"],
                    "protein": protein,
                    "equation": rxn["equation"],
                    "rhea_id": rxn["rhea_id"],
                }

                # Uniprot reg data (cached)
                entry["uniprot_reg_data"] = fetch_uniprot_reg_data(gene["accession"])

                # NCBI seq (cached)
                entry["reg_protein_seq"] = fetch_protein_fasta(gene["accession"])

                # Protein2chemicals (cached)
                # Try to gather info for all genes first, then process?
                protein_data = protein2chemicals(gene["accession"])
                if isinstance(protein_data, dict) and "catalysis" in protein_data:
                    ligand_names += protein_data["catalysis"].split(" ")

                # Filter ligands
                not_ligands = {"H2O", "+", "-", "=", "A", "AH2", "H(+)", "NADPH", "NADH", "NADP(+)", "NAD(+)", "2", "H(+)in", "H(+)out"}
                unique_ligands = [i for i in set(ligand_names) if i not in not_ligands]
                entry["alt_ligands"] = unique_ligands

                reg_data.append(entry)

    return reg_data
