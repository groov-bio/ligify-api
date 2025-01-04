from predict.chemical2enzymes import fetch_reactions, fetch_genes, filter_genes
from predict.enzymes2operons import pull_regulators
from predict.accID2operon import acc2OperonList
from predict.rank import calculate_rank
import requests
from pprint import pprint


def fetch_data(InChiKey, filters):
    metrics = {}

    # Enzyme reaction databases and their ligands:
        # RheaDB:   13,972
        # KEGG:     19,409
        # BRENDA:   173,436
        # MetaCyc:  19,400

    # TFBMiner gets reactions from KEGG via the REST API:
        # https://rest.kegg.jp/get/cpd:C00511
            # need to input KEGG compound ID, so we'd need to convert smiles into that
            # this will return EC numbers: 1.3.1.125       2.8.3.12        3.5.1.4         3.5.5.7  
            # from there you can use those EC numbers in new queries to get genes:
        # https://rest.kegg.jp/get/2.8.3.12
            # This will return the KEGG IDs of associated genes, if any:
            # EH206_20600, EH206_20605, Sant_2230, Sant_2231, ...
            # from here you can use each gene in a new query to get the ncbi ID:
        # https://rest.kegg.jp/get/xct:J151_00486
            # this returns a bunch of info, including the ncbi ID: AJD66958
            # at this point, the workflow merges with the parallel RHEA workflow


    # FETCH REACTIONS
    reactions = fetch_reactions(
        InChiKey=InChiKey, max_reactions=filters["max_reactions"]
    )
    total_rxns = len(reactions["rxn_data"])
    metrics["RHEA Reactions"] = total_rxns
    
    print(len(reactions["rxn_data"]) == 0)

    if total_rxns > 0:
        # FETCH GENES
        counter = 0
        for i in reactions["rxn_data"]:
            associated_proteins = fetch_genes(
                i["rhea_id"], filters["reviewed"], filters["proteins_per_reaction"]
            )
            i["proteins"] = associated_proteins
            counter += 1

        metrics["Total genes"] = sum(
            [len(i["proteins"]) for i in reactions["rxn_data"]]
        )

        # Filter homologous genes
        reactions = filter_genes(reactions, lineage_filter_name=filters["lineage"])
        metrics["Filtered genes"] = sum(
            [len(i["proteins"]) for i in reactions["rxn_data"]]
        )

        # FETCH OPERONS
        if len(reactions["rxn_data"]) == 0:
            print("No enzymes found")
            raise Exception(f"No enzymes found for {InChiKey}")
        else:
            # operon_counter = 0
            operon_list_entries = {}

            for rxn in reactions["rxn_data"]:
                if rxn["proteins"]:
                    for protein in rxn["proteins"]:
                        refseq_id = protein["enzyme"]["ncbi_id"]
                        if refseq_id is not None:
                            operon_list_entries[refseq_id] = None
                        # if (
                        #     refseq_id is not None
                        #     # and len(operon_list_entries.keys()) > 0
                        #     # <= filters["max_operons"]
                        # ):
                        #     operon_list_entries[refseq_id] = None

            acc2OperonListResult = acc2OperonList(operon_list_entries)

            metrics["Total operons"] = len(acc2OperonListResult.keys())

            # Attach result to context
            for rxn in reactions["rxn_data"]:
                if rxn["proteins"]:
                    for protein in rxn["proteins"]:
                        refseq_id = protein["enzyme"]["ncbi_id"]
                        protein["context"] = acc2OperonListResult[refseq_id]

            # FETCH REGULATORS

            if reactions is None:
                return None, None

            else:
                # This is where all of the display data is created
                counter = 0
                regulators = []

                for rxn in reactions["rxn_data"]:
                    for protein in rxn["proteins"]:
                        regs = pull_regulators(protein, rxn)
                        for r in regs:
                            regulators.append(r)
                        counter += 1

                metrics["Total regulators"] = len(regulators)

                # Filter out duplicate regulators
                refseq_ids = []
                filtered_regulators = []
                for i in regulators:
                    if i["refseq"] not in refseq_ids:
                        filtered_regulators.append(i)
                        refseq_ids.append(i["refseq"])

                # Filter out regulators without a predicted promoter
                filtered_regulators = [
                    i
                    for i in filtered_regulators
                    if i["protein"]["context"]["promoter"] is not None
                ]

                # Create a rank for each regulator
                for r in filtered_regulators:
                    rank = calculate_rank(r)
                    r["rank"] = rank

                if filtered_regulators is None or len(filtered_regulators) == 0:
                    raise Exception(f"No regulators found for {InChiKey}")
                else:
                    return filtered_regulators, metrics

    else:
        raise Exception(f"No reaction data found for {InChiKey}")


if __name__ == "__main__":

    def get_reactions(chem:str):
        url = "https://rest.kegg.jp/get/cpd:"
        response = requests.get(url+chem)

        if response.ok:
            try:
                for line in response.text.rstrip().split("\n"):
                    section = line[:12].strip()
                    if section == "ENZYME":
                        index = response.text.rstrip().split("\n").index(line)
                        reactions = line[12:].split(" ")
                        for line in response.text.rstrip().split("\n")[int(index)+1:]:
                            # If the current line is not part of the REACTION 
                            # section then the processing ends.
                            if line[:12].strip() != (""):
                                break
                            reactions_ = line[12:].split(" ")
                            reactions.extend(reactions_)
                # Eliminates erroneous IDs resulting from whitespaces.
                reactions = [reaction for reaction in reactions if reaction != ""]
                return reactions
            except UnboundLocalError:
                return []


    def get_genes(reaction):
        url = "https://rest.kegg.jp/get/"
        response = requests.get(url+reaction)        
        if response.ok:
            try:
                genes = []
                for line in response.text.rstrip().split("\n"):
                    section = line[:12].strip()
                    if section == "GENES":
                        index = response.text.rstrip().split("\n").index(line)
                        for line in response.text.rstrip().split("\n")[int(index)+1:]:
                            # If the current line is not part of the GENES 
                            # section then the processing ends.
                            if line[:12].strip() != (""):
                                break
                            genes_ = line[12:].split(" ")[0].lower() + line[12:].split(" ")[1]
                            #print(reactions_)
                            genes.append(genes_)
                # Eliminates erroneous IDs resulting from whitespaces.
                genes = [gene for gene in genes if gene != ""]
                return genes
            except UnboundLocalError:
                return []        


    def get_ncbiID(gene):
        url = "https://rest.kegg.jp/get/"
        response = requests.get(url+gene)        
        if response.ok:
            try:
                for line in response.text.rstrip().split("\n"):
                    section = line[:12].strip()
                    if section == "DBLINKS":
                        index = response.text.rstrip().split("\n").index(line)
                        lines = []
                        for line in response.text.rstrip().split("\n")[int(index):]:
                            # If the current line is not part of the DBLINKS 
                            # section then the processing ends.
                            if line[:7] == "DBLINKS":
                                lines.append(line)
                            elif line[:12].strip() != (""):
                                break
                            else: 
                                lines.append(line)
                            # process the DBLINKS-specific lines
                        for line in lines:
                            index = 0
                            found = False
                            for i in line.split(" "):
                                if i == 'NCBI-ProteinID:':
                                    found = True
                                    break
                                else:
                                    index += 1
                            if found:
                                return line.split(" ")[index+1]
            except UnboundLocalError:
                return []  


    reactions = get_reactions('C00511')
    genes = [k for i in reactions for k in get_genes(i)]
    print(len(genes))
    ncbiIDs = []
    for i in genes:
        IDs = get_ncbiID(i)
        if IDs is not None:
            for k in IDs:
                ncbiIDs.append(k)
    # ncbiIDs = [k for i in genes for k in get_ncbiID(i) if k is not None]
    print(ncbiIDs)
    #ncbiIDs = [get_ncbiID(i) for i in genes]
    #pprint(ncbiIDs)
    # genes = get_genes(reactions[0])
    # print(genes)
    # id = get_ncbiID(genes[0])
    # print(id)