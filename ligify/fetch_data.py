from predict.chemical2enzymes import fetch_reactions, fetch_genes, filter_genes
from predict.enzymes2operons import pull_regulators
from predict.accID2operon import acc2OperonList
from predict.rank import calculate_rank


def fetch_data(InChiKey, filters):
    metrics = {}

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
