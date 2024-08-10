from dotenv import load_dotenv

from ligify.fetch_data import fetch_data
from ligify.predict.pubchem import get_inchikey, get_name


if __name__ == "__main__":
    load_dotenv()
    smiles = "C=CC(=O)[O-]"
    chemical_name = get_name(smiles, "smiles")
    InChiKey = get_inchikey(smiles, "smiles")
    chemical = {"name": chemical_name, "smiles": smiles, "InChiKey": InChiKey}

    filters = {
        # Number
        "max_reactions": 20,
        # Number
        "proteins_per_reaction": 20,
        # bool
        "reviewed": True,
        # "Domain", "Phylum", "Class", "Order", "Family", "Genus", "None"
        "lineage": 'Family',
        # Number
        "max_operons": 20,
        # Number
        "max_alt_chems": 10
    }

    regulators, metrics = fetch_data(chemical['InChiKey'], filters)


    # {'RHEA Reactions': 5, 'Total genes': 8, 'Filtered genes': 4, 'Total operons': 4, 'Total regulators': 4}
    # print(regulators)
    print(metrics)