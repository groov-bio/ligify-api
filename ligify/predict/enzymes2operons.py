import json
import re
import requests


def protein2chemicals(accession: str):
    url = (
        "https://rest.uniprot.org/uniprotkb/search?query=" + accession + "&format=json"
    )

    response = requests.get(url)
    if response.ok:
        protein = json.loads(response.text)

        if len(protein["results"]) > 0:
            if "comments" in protein["results"][0]:
                protein = protein["results"][0]["comments"]

                protein_data = {}

                # look for description
                response_function = [
                    i["texts"][0]["value"]
                    for i in protein
                    if i["commentType"] == "FUNCTION"
                ]
                if len(response_function) != 0:
                    protein_data["function"] = response_function[0]

                    # look for catalytic activity
                catalysis = [
                    i["reaction"]["name"]
                    for i in protein
                    if i["commentType"] == "CATALYTIC ACTIVITY"
                ]
                if len(catalysis) != 0:
                    protein_data["catalysis"] = catalysis[0]

                    # look for ligands
                try:
                    reaction = [
                        i["reaction"]["reactionCrossReferences"]
                        for i in protein
                        if i["commentType"] == "CATALYTIC ACTIVITY"
                    ]
                    if len(reaction) != 0:
                        LIGANDS = [
                            i["id"] for i in reaction[0] if i["database"] == "ChEBI"
                        ]
                        if len(LIGANDS) != 0:
                            protein_data["ligands"] = LIGANDS
                except Exception:
                    pass

                    # look for induction
                induction = [
                    i["texts"][0]["value"]
                    for i in protein
                    if i["commentType"] == "INDUCTION"
                ]
                if len(induction) != 0:
                    protein_data["induction"] = induction[0]

                    # look for induction
                pathway = [
                    i["texts"][0]["value"]
                    for i in protein
                    if i["commentType"] == "PATHWAY"
                ]
                if len(pathway) != 0:
                    protein_data["pathway"] = pathway[0]

                # add something to append all this metadata

                return protein_data
    else:
        response.raise_for_status()


def fetch_reg_protein_seq(accession: str):
    URL = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="
        + accession
        + "&rettype=fasta"
    )
    response = requests.get(URL)
    if response.ok:
        seq = "".join(i for i in response.text.split("\n")[1:])
        return seq
    else:
        print("Bad eFetch request " + str(response.status_code))
        return None


def fetch_uniprot_reg_data(accession: str):
    url = (
        "https://rest.uniprot.org/uniprotkb/search?query=" + accession + "&format=json"
    )
    response = requests.get(url)

    try:
        data = json.loads(response.text)["results"][0]

        dois = []
        for j in data["references"]:
            if "citationCrossReferences" in j["citation"]:
                for k in j["citation"]["citationCrossReferences"]:
                    if k["database"] == "DOI":
                        doi = k["id"]
                        title = j["citation"]["title"]
                        break
                dois.append({"doi": doi, "title": title})

        regulator = {
            "annotation": data["features"][0]["description"],
            "id": data["primaryAccession"],
            "references": dois,
            "length": data["sequence"]["length"],
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
    regulator = re.compile(r"regulator|repressor|activator")

    reg_data = []

    ligand_names = []
    if "context" in protein.keys():
        if protein["context"] != "EMPTY":
            operon = protein["context"]["operon"]
            for gene in operon:
                if "description" in gene.keys():
                    if regulator.search(gene["description"]):
                        entry = {
                            "refseq": gene["accession"],
                            "annotation": gene["description"],
                            "protein": protein,
                            "equation": rxn["equation"],
                            "rhea_id": rxn["rhea_id"],
                        }

                        ### This is where a Uniprot API query goes to fetch more info on the regulator.
                        entry["uniprot_reg_data"] = fetch_uniprot_reg_data(
                            gene["accession"]
                        )

                        ### NCBI is queried for more info on the regulator.
                        # This is a fail-safe, since sometimes no Uniprot ID is associated with the RefSeq identifier
                        entry["reg_protein_seq"] = fetch_reg_protein_seq(
                            gene["accession"]
                        )

                        # Fetch possible alternative inducer molecules associated with the operon





                        # Fetch possible alternative inducer molecules associated with the operon
                        ligand_ids = []
                        for gene in operon:
                            protein_data = protein2chemicals(gene["accession"])
                            if isinstance(protein_data, dict):
                                if "ligands" in protein_data.keys():
                                    for l in protein_data['ligands']:
                                        ligand_ids.append(l)
                        unique_ligand_ids = list(set(ligand_ids))
                        # Blacklisted ligands
                        not_ligand_ids = [
                            "CHEBI:15378",    # H(+)
                            "CHEBI:15377",    # H2O
                            "CHEBI:16474",    # NADPH
                            "CHEBI:16908",    # NADH
                            "CHEBI:57945",    # NADH(-2)
                            "CHEBI:57540",    # NAD(-1)
                            "CHEBI:18009",    # NADP(+)
                            "CHEBI:15846",    # NAD(+)
                        ]
                        unique_ligand_ids = [
                            i for i in unique_ligand_ids if i not in not_ligand_ids
                        ]

                        # function to get name and smiles
                        def get_smiles_and_name(input):
                            URL = (
                                "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
                                + str(input)
                                + "/property/IsomericSMILES,IUPACName/JSON"
                            )
                            response = requests.get(URL)
                            if response.ok:
                                data = json.loads(response.text)
                                name = data["PropertyTable"]["Properties"][0]["IUPACName"]
                                smiles = data["PropertyTable"]["Properties"][0]["IsomericSMILES"]

                                return {"name":name, "smiles":smiles}

                        ligands = []
                        for i in unique_ligand_ids:
                            ligands.append(get_smiles_and_name(i))


                        entry["candidate_ligands"] = ligands









                        for gene in operon:
                            protein_data = protein2chemicals(gene["accession"])
                            if isinstance(protein_data, dict):
                                if "catalysis" in protein_data.keys():
                                    ligand_names += protein_data["catalysis"].split(" ")
                        unique_ligands = list(set(ligand_names))
                        # Blacklisted ligands
                        not_ligands = [
                            "H2O",
                            "+",
                            "-",
                            "=",
                            "A",
                            "AH2",
                            "H(+)",
                            "NADPH",
                            "NADH",
                            "NADP(+)",
                            "NAD(+)",
                            "2",
                            "H(+)in",
                            "H(+)out",
                        ]
                        unique_ligands = [
                            i for i in unique_ligands if i not in not_ligands
                        ]

                        entry["alt_ligands"] = unique_ligands

                        reg_data.append(entry)

    return reg_data


if __name__ == "__main__":
    

    # function to get name and smiles
    def get_smiles_and_name(input):
        URL = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            + str(input)
            + "/property/IsomericSMILES,IUPACName/JSON"
        )
        response = requests.get(URL)
        if response.ok:
            data = json.loads(response.text)
            name = data["PropertyTable"]["Properties"][0]["IUPACName"]
            smiles = data["PropertyTable"]["Properties"][0]["IsomericSMILES"]

            return {"name":name, "smiles":smiles}


    print(get_smiles_and_name("CHEBI:15378"))