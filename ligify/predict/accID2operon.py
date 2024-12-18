import json
import os
import requests
import re
import xmltodict

# TODO:
# Return a legit error message for the frontend if an error comes up

# OLD VERSION

# def NC2genome(genome_id, startPos, stopPos):

#     base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore"
#     response = requests.get(base_url+"&id="+str(genome_id)+"&seq_start="+str(startPos)+"&seq_stop="+str(stopPos)+"&rettype=fasta_cds_aa")
#     # old script that fetches the entire genome
#     # response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='+str(genome_id)+'&rettype=fasta_cds_aa')

#     if response.ok:
#         genome = response.text.split("\n")
#         return genome


def NC2genome(genome_id, operon):
    base_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&api_key={os.getenv('NcbiApiKey')}"
    startPos = operon[0]["start"]
    stopPos = operon[-1]["stop"]
    response = requests.get(
        base_url
        + "&id="
        + str(genome_id)
        + "&seq_start="
        + str(startPos)
        + "&seq_stop="
        + str(stopPos)
        + "&rettype=fasta"
    )

    if response.ok:
        genome = response.text.split("\n")
        genome = "".join(i for i in genome[0:-2] if i[0] != ">")

        ### GENOME FRAGMENT ANNOTATION FUNCTION ###

        ### This deals with one-sided gene overlaps (beginning or end)

        ### It does NOT YET deal with double-sided overlaps (begining AND end)

        out = {}
        counter = 0
        for index in range(0, len(operon)):
            # reset overlap seq
            overlap_seq = ""

            # if you're not at the end...
            if index != len(operon) - 1:
                # if END of gene overlaps with START of next gene...
                if operon[index + 1]["start"] < operon[index]["stop"]:
                    # truncated gene
                    gene_seq = genome[
                        operon[index]["start"] - startPos : operon[index + 1]["start"]
                        - startPos
                    ]
                    # overlap region
                    overlap_seq = genome[
                        operon[index + 1]["start"] - startPos : operon[index]["stop"]
                        - startPos
                        + 1
                    ]

                # if you're not at the beginning...
                elif index != 0:
                    # if START of gene overlaps with END of prior gene...
                    if operon[index - 1]["stop"] > operon[index]["start"]:
                        # truncated gene
                        gene_seq = genome[
                            operon[index - 1]["stop"] - startPos + 1 : operon[index][
                                "stop"
                            ]
                            - startPos
                            + 1
                        ]
                    else:
                        # full gene
                        gene_seq = genome[
                            operon[index]["start"] - startPos : operon[index]["stop"]
                            - startPos
                            + 1
                        ]

                # if you're at the beginning
                elif index == 0:
                    # full gene
                    gene_seq = genome[
                        operon[index]["start"] - startPos : operon[index]["stop"]
                        - startPos
                        + 1
                    ]

            # if you ARE at the end...
            elif index == len(operon) - 1:
                # see if START of gene overlaps with END of prior gene
                if operon[index - 1]["stop"] > operon[index]["start"]:
                    # truncated gene
                    gene_seq = genome[
                        operon[index - 1]["stop"] - startPos + 1 : operon[index]["stop"]
                        - startPos
                        + 1
                    ]
                else:
                    # full gene
                    gene_seq = genome[operon[index]["start"] - startPos :]

            # Append the gene sequence
            if str(operon[index]["direction"]) == "+":
                out["gene" + str(counter) + "fwd"] = gene_seq
            else:
                out["gene" + str(counter)] = gene_seq

            # Append the overlap sequence
            if len(overlap_seq) > 0:
                out["overlap" + str(counter)] = overlap_seq

            # Append the spacer
            try:
                spacer_seq = genome[
                    operon[index]["stop"] - startPos + 1 : operon[index + 1]["start"]
                    - startPos
                ]
            except Exception:
                spacer_seq = genome[operon[index]["stop"] - startPos + 1 : len(genome)]
            if len(spacer_seq) != 0:
                out["spacer" + str(counter)] = spacer_seq

            counter += 1

        reconstruct = "".join(i for i in out.values())
        if reconstruct == genome:
            genome_reassembly_match = True
        else:
            genome_reassembly_match = False

        return out, genome_reassembly_match


def getGenes(genome_id, startPos, stopPos):
    # Fetch the genome fragment
    base_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&api_key={os.getenv('NcbiApiKey')}"
    try:
        response = requests.get(
            base_url
            + "&id="
            + str(genome_id)
            + "&seq_start="
            + str(startPos - 10000)
            + "&seq_stop="
            + str(stopPos + 10000)
            + "&rettype=fasta_cds_aa"
        )
        genome = response.text.split("\n")
    except Exception:
        try:
            response = requests.get(
                base_url
                + "&id="
                + str(genome_id)
                + "&seq_start="
                + str(startPos - 5000)
                + "&seq_stop="
                + str(stopPos + 5000)
                + "&rettype=fasta_cds_aa"
            )
            genome = response.text.split("\n")
        except Exception:
            try:
                response = requests.get(
                    base_url
                    + "&id="
                    + str(genome_id)
                    + "&seq_start="
                    + str(startPos)
                    + "&seq_stop="
                    + str(stopPos + 5000)
                    + "&rettype=fasta_cds_aa"
                )
                genome = response.text.split("\n")
            except Exception:
                try:
                    response = requests.get(
                        base_url
                        + "&id="
                        + str(genome_id)
                        + "&seq_start="
                        + str(startPos - 5000)
                        + "&seq_stop="
                        + str(stopPos)
                        + "&rettype=fasta_cds_aa"
                    )
                    genome = response.text.split("\n")
                except Exception:
                    print("error fetching the genome fragment")
                    return None, None

    re1 = re.compile(str(startPos))
    re2 = re.compile(str(stopPos))
    geneIndex = 0
    regIndex = None
    genes = []
    for i in genome:
        if len(i) != 0:
            if i[0] == ">":
                if re1.search(i):
                    if re2.search(i):
                        regIndex = geneIndex
                geneIndex += 1
                genes.append(i)
    if regIndex is None:
        print("regulator not found in genome")
        return None, None
    else:
        return genes, regIndex


def fasta2MetaData(fasta):
    metaData = {}
    regulator = fasta.split(" [")

    for i in regulator:
        if i[:10] == "locus_tag=":
            metaData["alias"] = i[10:-1]
        elif i[:8] == "protein=":
            metaData["description"] = i[8:-1].replace("'", "")
        elif i[:11] == "protein_id=":
            metaData["accession"] = i[11:-1]
        elif i[:9] == "location=":
            if i[9:20] == "complement(":
                metaData["direction"] = "-"
                location = i[20:-2]
                location = location.split("..")
                metaData["start"] = int(re.sub("\D", "", location[0]))
                metaData["stop"] = int(re.sub("\D", "", location[1]))
            else:
                metaData["direction"] = "+"
                location = i[9:-1]
                location = location.split("..")
                metaData["start"] = int(re.sub("\D", "", location[0]))
                metaData["stop"] = int(re.sub("\D", "", location[1]))

    if "accession" not in metaData.keys():
        metaData["accession"] = ""

    return metaData


def getOperon(allGenes, index, seq_start, strand):
    """
    Rules for inclusion/exclusion of genes from operon:
        - always take immediately adjacent genes
        - if query gene is in same direction as regulator, include it.
        - if query gene is expressed divergently from regulator,
                grab all adjacent genes that are expressed divergently (change strand direction for next genes)
        - if query gene is expressed divergently from a co-transcribed neighbor of the regulaor,
                grab that gene. (it may be another regulator. Important to know).
        - if query gene direction converges with regulator, exclude it.
    """

    def getGene(geneStrand, direction, nextGene, geneList, index):
        while geneStrand == nextGene["direction"]:
            if direction == "+":
                nextIndex = index + 1
            elif direction == "-":
                nextIndex = index - 1

            try:
                nextGene = fasta2MetaData(allGenes[nextIndex])

                if (
                    abs(seq_start - nextGene["start"]) > 8000
                ):  # added this. break if too far away
                    break

                elif (
                    geneStrand == "-"
                    and nextGene["direction"] == "+"
                    and direction == "+"
                ):
                    geneList.append(nextGene)
                elif (
                    geneStrand == "+"
                    and nextGene["direction"] == "-"
                    and direction == "-"
                ):
                    geneList.append(nextGene)
                elif geneStrand == nextGene["direction"]:
                    geneList.append(nextGene)
                index = nextIndex
            except Exception:
                break

    geneStrand = strand

    # attempt to get downstream genes, if there are any genes downstream
    try:
        indexDOWN = index - 1
        downGene = fasta2MetaData(allGenes[indexDOWN])
        # if seq_start > downGene['start']:
        if strand == "+" and downGene["direction"] == "-":
            geneStrand = downGene["direction"]

        downgenes = [downGene]
        getGene(geneStrand, "-", downGene, downgenes, indexDOWN)

        geneArray = list(reversed(downgenes))
    except Exception:
        geneArray = []

    geneArray.append(fasta2MetaData(allGenes[index]))
    regulatorIndex = len(geneArray) - 1

    geneStrand = strand

    # attempt to get upstream genes, if there are any genes upstream
    try:
        indexUP = index + 1
        upGene = fasta2MetaData(allGenes[indexUP])
        # if seq_start > upGene['start']:
        if strand == "-" and upGene["direction"] == "+":
            geneStrand = upGene["direction"]

        geneArray.append(upGene)

        getGene(geneStrand, "+", upGene, geneArray, indexUP)
    except Exception:
        return geneArray, regulatorIndex

    return geneArray, regulatorIndex


def predict_promoter(operon, regIndex, genome_id):
    if operon[regIndex]["direction"] == "+":
        queryGenes = list(reversed(operon[0:regIndex]))
        index = regIndex
        if len(queryGenes) == 0:
            # print("WARNING: Tiny operon with too few genes. This entry will be omitted.")
            return
        for i in queryGenes:
            if i["direction"] == "-":
                startPos = i["stop"]
                stopPos = operon[index]["start"]
                regType = 1
                break
            else:
                start = operon[regIndex - 1]["stop"]
                stop = operon[regIndex]["start"]
                testLength = int(stop) - int(start)
                # Setting this to 100bp is somewhat arbitrary.
                # Most intergenic regions >= 100bp. May need to tweak.
                if testLength > 100:
                    startPos = start
                    stopPos = stop
                    regType = 2
                    break
                else:
                    if index == 1:
                        # print('WARNING: Reached end of operon. This entry will be omitted')
                        return None
                    index -= 1

    elif operon[regIndex]["direction"] == "-":
        queryGenes = operon[regIndex + 1 :]
        index = regIndex
        if len(queryGenes) == 0:
            # print("WARNING: Tiny operon with too few genes. This entry will be omitted.")
            return
        for i in queryGenes:
            if i["direction"] == "+":
                stopPos = i["start"]
                startPos = operon[index]["stop"]
                regType = 1
                break
            else:
                # Counterintunitive use of "stop"/"start" ...
                # Start < stop always true, regardless of direction
                start = operon[regIndex]["stop"]
                stop = operon[regIndex + 1]["start"]
                testLength = int(stop) - int(start)
                if testLength > 100:
                    startPos = start
                    stopPos = stop
                    regType = 2
                    break
                else:
                    if index == len(operon) - 2:
                        # print('WARNING: Reached end of operon. This entry will be omitted')
                        return None
                    else:
                        index += 1

    URL = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="
        + str(genome_id)
        + "&seq_start="
        + str(startPos)
        + "&seq_stop="
        + str(stopPos)
        + "&strand=1&rettype=fasta"
        + f"&api_key={os.getenv('NcbiApiKey')}"
    )
    response = requests.get(URL)

    if response.ok:
        intergenic = response.text
        output = ""
        for i in intergenic.split("\n")[1:]:
            output += i
        if len(output) <= 1000:
            return {"regulated_seq": output[1:-1], "reg_type": regType}
        else:
            # TODO: This is a potential failure mode!!!
            # print('WARNING: Intergenic region is over 800bp')
            return None
    else:
        print(f"Status Code: {response.status_code}")
        print(f"Reason: {response.reason}")
        print(f"Response Text: {response.text}")
        print(f"Response Headers: {response.headers}")
        print("FATAL: Bad eFetch request")
        return None

        # 800bp cutoff for an inter-operon region.
        # A region too long makes analysis fuzzy and less accurate.


def acc2MetaDataList(access_ids):
    base_url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=ipg&api_key={os.getenv('NcbiApiKey')}"
    )
    id_list = list(access_ids.keys())
    for id in id_list:
        access_ids[id] = None  # Initialize to None
        base_url += f"&id={id}"

    result = requests.get(base_url)
    if result.status_code != 200:
        print("non-200 HTTP response. eFetch failed")
        print(f"Status Code: {result.status_code}")
        print(f"Reason: {result.reason}")
        print(f"Response Text: {result.text}")
        print(f"Response Headers: {result.headers}")
        print("non-200 HTTP response. eFetch failed")
        return access_ids  # Return early if the request fails

    parsed = xmltodict.parse(result.text)

    if "IPGReport" in parsed["IPGReportSet"].keys():
        # Ensure IPGReport is a list
        ipg_reports = parsed["IPGReportSet"]["IPGReport"]
        if not isinstance(ipg_reports, list):
            ipg_reports = [ipg_reports]

        ids_set_to_none = set()

        for entry in ipg_reports:
            input_acc = entry.get("@product_acc")
            if not input_acc:
                # Cannot tie back to access_ids
                print("No '@product_acc' in entry, cannot tie back to access_ids")
                continue

            id = input_acc

            if "ProteinList" in entry:
                protein = entry["ProteinList"]["Protein"]
                if isinstance(protein, list):
                    protein = protein[0]

                if "CDSList" not in protein:
                    access_ids[id] = None  # Explicitly set to None
                    ids_set_to_none.add(id)
                    continue

                CDS = protein["CDSList"]["CDS"]
                if isinstance(CDS, list):
                    CDS = CDS[0]

                proteinDict = {
                    "accver": CDS["@accver"],
                    "start": CDS["@start"],
                    "stop": CDS["@stop"],
                    "strand": CDS["@strand"],
                }
                access_ids[id] = proteinDict
            else:
                access_ids[id] = None  # Explicitly set to None
                ids_set_to_none.add(id)
    else:
        # Handle error if IPGReport is missing
        raise Exception("No IPGReport in response")

    # Fallback for IDs not explicitly set to None
    for id in access_ids:
        if access_ids[id] is None and id not in ids_set_to_none:
            # Use the fallback function
            access_ids[id] = acc2MetaData(id)

    return access_ids


def acc2MetaData(access_id: str):
    result = requests.get(
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={access_id}&rettype=ipg&api_key={os.getenv('NcbiApiKey')}"
    )
    if result.status_code != 200:
        print(f"Status Code: {result.status_code}")
        print(f"Reason: {result.reason}")
        print(f"Response Text: {result.text}")
        print(f"Response Headers: {result.headers}")
        print("non-200 HTTP response. eFetch failed")
        return None  # Return None on failure

    parsed = xmltodict.parse(result.text)

    if "IPGReport" in parsed["IPGReportSet"].keys():
        if "ProteinList" in parsed["IPGReportSet"]["IPGReport"]:
            protein = parsed["IPGReportSet"]["IPGReport"]["ProteinList"]["Protein"]

            if isinstance(protein, list):
                protein = protein[0]

            if "CDSList" not in protein:
                return None  # Return None if no CDSList

            CDS = protein["CDSList"]["CDS"]

            # CDS is a list if there is more than 1 CDS returned, otherwise it's a dictionary
            if isinstance(CDS, list):
                CDS = CDS[0]

            proteinDict = {
                "accver": CDS["@accver"],
                "start": CDS["@start"],
                "stop": CDS["@stop"],
                "strand": CDS["@strand"],
            }

            return proteinDict
        else:
            return None
    else:
        return None



def acc2OperonList(operon_list_entries):
    metaData = acc2MetaDataList(operon_list_entries)
    operon_dict = {}

    for key, value in metaData.items():
        if value is not None:
            genes, index = getGenes(
                value["accver"], int(value["start"]), int(value["stop"])
            )

            if index is not None:
                enzyme = fasta2MetaData(genes[index])

                operon, regIndex = getOperon(
                    genes, index, enzyme["start"], enzyme["direction"]
                )
                operon_sequence, reassembly_match = NC2genome(value["accver"], operon)
                promoter = predict_promoter(operon, regIndex, value["accver"])

                data = {
                    "operon": operon,
                    "enzyme_index": regIndex,
                    "enzyme_direction": enzyme["direction"],
                    "operon_seq": operon_sequence,
                    "promoter": promoter,
                    "reassembly_match": reassembly_match,
                    "genome": value["accver"],
                }

                # OLD
                # data = {"operon": operon, "enzyme_index": regIndex, "genome": metaData["accver"] }

                operon_dict[key] = data
            else:
                operon_dict[key] = "EMPTY"
        else:
            operon_dict[key] = "EMPTY"

    return operon_dict


def acc2operon(accession):
    metaData = acc2MetaData(accession)

    if metaData != "EMPTY":
        genes, index = getGenes(
            metaData["accver"], int(metaData["start"]), int(metaData["stop"])
        )

        if index is not None:
            enzyme = fasta2MetaData(genes[index])

            operon, regIndex = getOperon(
                genes, index, enzyme["start"], enzyme["direction"]
            )
            operon_sequence, reassembly_match = NC2genome(metaData["accver"], operon)
            promoter = predict_promoter(operon, regIndex, metaData["accver"])

            data = {
                "operon": operon,
                "enzyme_index": regIndex,
                "enzyme_direction": enzyme["direction"],
                "operon_seq": operon_sequence,
                "promoter": promoter,
                "reassembly_match": reassembly_match,
                "genome": metaData["accver"],
            }

            # OLD
            # data = {"operon": operon, "enzyme_index": regIndex, "genome": metaData["accver"] }

            return data
        else:
            return "EMPTY"
    else:
        return "EMPTY"


if __name__ == "__main__":
    # metaData = acc2MetaData('WP_011046214.1')

    # if metaData != "EMPTY":
    #     genes, index = getGenes(metaData["accver"], int(metaData["start"]), int(metaData["stop"]))
    #     print(genes)
    #     print(index)

    # print(acc2MetaData('NP_414848.1'))

    print(acc2MetaDataList({'NP_414848.1': None, 'NP_390983.2': None, 'WP_011369019.1': None, 'WP_013240889.1': None, 'WP_041640381.1': None}))

    # print(acc2operon('WP_011046214.1'))
    # acc2operonList(['WP_011046214.1', 'ADT64689.1', 'WP_011336734.1','WP_069687654.1'])
    # acc2operonList(['WP_011046214.1', 'WP_011336734.1'])
    # pprint(acc2operonList(['WP_069687654.1', 'WP_011046214.1', 'ADT64689.1', 'WP_011336734.1']))
