from dnachisel import (
    CodonOptimize,
    reverse_translate,
    DnaOptimizationProblem,
    AvoidPattern,
    EnforceGCContent,
    EnforceTranslation,
)


def codon_opt(protein_seq: str):
    # Create a random DNA seq given the protein seq. Append a stop codon.
    try:
        protein_dna_seq = reverse_translate(protein_seq + "*")
    except Exception as e:
        print("reverse_translate exception")
        print(e)

    # DEFINE THE OPTIMIZATION PROBLEM
    try:
        problem = DnaOptimizationProblem(
            sequence=protein_dna_seq,
            constraints=[
                AvoidPattern("BsaI_site"),
                EnforceGCContent(mini=0.35, maxi=0.65, window=50),
                EnforceTranslation(location=(0, len(protein_dna_seq))),
            ],
            objectives=[
                CodonOptimize(species="e_coli", location=(0, len(protein_dna_seq)))
            ],
        )
    except Exception as e:
        print("DnaOptimizationProblem exception")
        print(e)

    # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE

    try:
        problem.resolve_constraints()
    except Exception as e:
        print("resolve_constraints exception")
        print(e)
    try:
        problem.optimize()
    except Exception as e:
        print("optimize exception")
        print(e)

    # GET THE FINAL SEQUENCE (AS STRING OR ANNOTATED BIOPYTHON RECORDS)

    final_sequence = problem.sequence  # string
    return final_sequence


if __name__ == "__main__":
    codon_opt(
        "MTTIRWRRMSIHSERITLADSPLHWAHTLNGSMRTHFEVQRLERGRGAYLARSRFGAGELYSAIAPSQVLRHFNDQRNANEAEHSYLIQIRSGALGVASGGRKVILANGDCSIVDSRQDFTLSSNSSTQGVVIRFPVSWLGAWVSNPEDLIARRVDAEIGWGRALSASVSNLDPLRIDDLGSNVNSIAEHVAMLISLASSAVSSEDGGVALRKMREVKRVLEQSFADANLEPESVSSQLGISKRYLHYVFAACGTTFGRELLEIRLGKAYRMLCATSGSGAVLKVAMSSGFSDSSHFSKKFKERYGVSPVSLVRQA"
    )
