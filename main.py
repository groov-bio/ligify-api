from pprint import pprint
from dotenv import load_dotenv
import json

from ligify.fetch_data import fetch_data
from ligify.genbank.create_genbank import create_genbank
from ligify.predict.pubchem import get_inchikey, get_name

# TODO - accept smiles or inchikey prop

def create_plasmid(regulators, chemical):
    for regulator in regulators:
        result = create_genbank(regulator['refseq'], chemical, regulator['protein']['context']['promoter']['regulated_seq'], regulator['reg_protein_seq'])
        regulator['plasmid_sequence'] = str(result)
    
    return regulators

if __name__ == "__main__":
    load_dotenv()
    smiles = "C=CC(=O)[O-]"
    input_type = "smiles"
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
    try:
        chemical_name = get_name(smiles, input_type)
        InChiKey = get_inchikey(smiles, input_type)
        print(chemical_name)
        chemical = {"name": chemical_name, "smiles": smiles, "InChiKey": InChiKey}

        regulators, metrics = fetch_data(chemical['InChiKey'], filters)
        regulators = create_plasmid(regulators, chemical_name)
        # for regulator in regulators:
        #     print(regulator['refseq'])
        with open('data.json', 'w', encoding='utf-8') as f:
            json.dump(regulators, f, ensure_ascii=False, indent=4)
        # pprint(regulators)
        # format_results(data_column, chemical["name"])
        # {'RHEA Reactions': 5, 'Total genes': 8, 'Filtered genes': 4, 'Total operons': 4, 'Total regulators': 4}
        print(metrics)
    except Exception as e:
        print(e)

# TODO 
# Hookup create_genbank with the data returned from fetch_data
# Walk through data.json and see if anything is no longer relevant
# Cleanse create_genbank of streamlit
