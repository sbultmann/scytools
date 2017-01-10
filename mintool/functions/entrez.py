def query_ncbi(gene_name, org_name):
    """ Takes Gene Name and Organism to get related genes and isoform
            information @ Entrez """
    from Bio import Entrez
    import json
    import re
    Entrez.email = "bultmann@bio.lmu.de"
    Entrez.tool = "mintool"

    handle = Entrez.esearch(db="gene", term="%s[Gene Name] AND %s[Organism]" % (gene_name,org_name))
    record = Entrez.read(handle)
    ids = record["IdList"][0]

    handle = Entrez.efetch(db="gene", id=ids, rettype="gb", retmode="xml")
    records = Entrez.parse(handle)
    result = []
    for record in records:
        x = record['Entrezgene_locus'][0]['Gene-commentary_products']
        gene_ref = record['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']
        gene_syn = record['Entrezgene_gene']['Gene-ref']['Gene-ref_syn']
        locus_accession = record['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']\
            ['Seq-interval']['Seq-interval_id']['Seq-id']['Seq-id_gi']
        locus_toT =record['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']\
            ['Seq-interval']['Seq-interval_to']
        locus_fromT = record['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int'] \
            ['Seq-interval']['Seq-interval_from']
        locus_to = max([int(locus_toT), int(locus_fromT)])
        locus_from = min([int(locus_toT), int(locus_fromT)])

        locus_strand = str(record['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int'] \
            ['Seq-interval']['Seq-interval_strand']).split(":")[2].split("}")[0][2:-1]
        if locus_strand == "minus":
            locus_strand = 2
        else:
            locus_strand = 1
        locus_chr = record['Entrezgene_locus'][0]['Gene-commentary_label']
        try:
            locus_des = record['Entrezgene_summary']
        except:
            locus_des = "No description available"
        iso_info_list = []
        pxwidth = 300
        maxlenght = max([int(locus_to), int(locus_from)])
        width = max([int(locus_to), int(locus_from)]) - min([int(locus_to), int(locus_from)])
        ratio = pxwidth/width
        for i in x:
            if 'Gene-commentary_products' in i.keys() and "X" not in i['Gene-commentary_label']:
                exons = []
                real_exons=[]
                isoform_name = i['Gene-commentary_label']
                isoform_accessionP = i['Gene-commentary_products'][0]['Gene-commentary_accession']
                isoform_accessionN = i['Gene-commentary_seqs'][0]['Seq-loc_whole']['Seq-id']['Seq-id_gi']
                isoform_locs = i['Gene-commentary_products'][0]['Gene-commentary_genomic-coords']\
                    [0]['Seq-loc_mix']['Seq-loc-mix']
                for loc in isoform_locs:
                    loc_strand = str(loc['Seq-loc_int']['Seq-interval']['Seq-interval_strand']).split(":")[2].split("}")[0][2:-1]
                    if loc_strand == "plus":
                        cor = 300
                    else:
                        cor = 0
                    loc_to = int(loc['Seq-loc_int']['Seq-interval']['Seq-interval_to'])
                    loc_from = int(loc['Seq-loc_int']['Seq-interval']['Seq-interval_from'])
                    real_exons.append(loc_to)
                    real_exons.append(loc_from)
                    loc_to = abs(int((float(loc['Seq-loc_int']['Seq-interval']['Seq-interval_to'])-maxlenght)*ratio)+cor)
                    loc_from = abs(int((float(loc['Seq-loc_int']['Seq-interval']['Seq-interval_from'])-maxlenght)*ratio)+cor)

                    #if loc_to > loc_from:
                    exon = [loc_to, loc_from]
                    #else:
                     #   exon = [loc_from, loc_to]
                    exons.append(exon)
                seqDict = extract_sequences(locus_accession,min(real_exons), max(real_exons), locus_strand)
                iso_info_list.append(json.dumps({"Name": isoform_name, "Accession_protein": isoform_accessionP,
                                                 "Accession_nucleotide": isoform_accessionN, "Exons": exons,
                                                 "Nterm": seqDict["Nterm"], "Cterm" : seqDict["Cterm"]}))

    result.append({"Gene_Name": gene_ref, "chr": re.findall(r"[0-9XY]+", locus_chr)[0],"Gene_Syn": gene_syn, "Accession": locus_accession, "To": locus_to,
                   "From": locus_from, "Strand": locus_strand, "Chr": locus_chr, "Description": locus_des,
                   "Isoforms": iso_info_list})

    return result


def extract_sequences(nuc_id, iso_start, iso_stop, strand):
    from Bio import Entrez
    from Bio import SeqIO
    Entrez.email = "bultmann@bio.lmu.de"
    Entrez.tool = "mintool"
    handle = Entrez.efetch(db="nucleotide", id=nuc_id, rettype="fasta", retmode="text",
                           seq_start=iso_start-499, seq_stop=iso_stop+501, strand=strand)
    record = SeqIO.read(handle, "fasta")
    Nterm = str(record.seq[0:1003])
    Cterm = str(record.seq[-1003:])

    return {"Nterm": Nterm, "Cterm" : Cterm}


def generate_grna(sequence, terminus):
    import re
    from Bio.Seq import Seq
    minTag = "ggtttgtctggtcaaccaccgcggtctcagtggtgtacggtacaaacc"
    if terminus == "seqnterm":
        modified_sequence = sequence[0:503] + minTag + sequence[503:]
    if terminus == "seqcterm":
        modified_sequence = sequence[0:500] + minTag + sequence[500:]

    positions = {}
    matches_fwd = re.finditer(r'(?=([ACGTacgt]{21}GG))', sequence)
    matches_rev = re.finditer(r'(?=(CC[ACGTacgt]{21}))', sequence)
    n = 1
    for match in matches_fwd:
        if match.start(1) > 450 and match.end(1) < 553:
            if match.group(1) not in modified_sequence:
                rating = -1
            else:
                rating = abs(match.end() - 500)
            positions["gRNA_" + str(n)] = {"start": match.start(1), "end": match.end(1), "sequence": match.group(1),
                                           "rating": rating, "strand": 1,
                                           "toligo": str(Seq(modified_sequence[427:627]).reverse_complement())}
            n += 1
    for match in matches_rev:
        if match.start(1) > 450 and match.end(1) < 553:
            if match.group(1) not in modified_sequence:
                rating = -1
            else:
                rating = abs(match.start() - 500)
            positions["gRNA_" + str(n)] = {"start": match.start(1), "end": match.end(1),
                                           "sequence": str(Seq(match.group(1)).reverse_complement()),
                                           "rating": rating, "strand": 2,
                                           "toligo": modified_sequence[427:627]}
            n += 1

    return positions


def design_primers(sequence, included_region, size_range):
    from primer3 import bindings
    primers = bindings.designPrimers(
        {
            'SEQUENCE_ID': 'MH1000',
            'SEQUENCE_TEMPLATE': sequence,
            'SEQUENCE_INCLUDED_REGION': included_region,
            'SEQUENCE_EXCLUDED_REGION': [400, 200]
        },
        {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE': size_range,
            'PRIMER_NUM_RETURN': 1
        })
    """parts needed:PRIMER_RIGHT_0_SEQUENCE, PRIMER_RIGHT_0, PRIMER_RIGHT_0_TM
                    PRIMER_LEFT_0_SEQUENCE, PRIMER_LEFT_0, PRIMER_LEFT_0_TM"""

    return {"left_pos": list(primers["PRIMER_LEFT_0"]), "left_tm": primers["PRIMER_LEFT_0_TM"],
            "left_seq": primers["PRIMER_LEFT_0_SEQUENCE"],
            "right_pos": list(primers["PRIMER_RIGHT_0"]), "right_tm": primers["PRIMER_RIGHT_0_TM"],
            "right_seq": primers["PRIMER_RIGHT_0_SEQUENCE"],
            "product_size": list(primers["PRIMER_RIGHT_0"])[0] - primers["PRIMER_LEFT_0"][0]
            }


def gRNA_oligos(seq):
    from Bio.Seq import Seq
    seq = Seq(seq)
    f_oligo = "caccg"+str(seq)[:-3]
    r_oligo = "aaac"+str(seq.reverse_complement())[3:]+"c"
    return (f_oligo, r_oligo)

#print(query_ncbi("Dnmt3a", "Mouse"))

#print(extract_sequences(372099101,20907520,20941453,2))
#print(extract_sequences(372099098,3835197,3907746,1))

#20908751 20908471 2 372099101
#print(design_primers("CTGAGCTCACCAGGTCAGCCCAGAGTCCCACTCTTTGAGAGCCCCTCTCTTCTCAGAGGGGCATTTATGGATGAATCAGATTGAATTCAGACTCAAGCTGCCTGCACCATCCAGTTCTCATGTGCAGGGCGCCTGTAGAGCACATCCAAGAGGCATTTCAAAGGTCTGGGGGTTTTGCTTTGTTTTGTTTTTTTTTACTGCCTTTGTTTCTGGGGAGTGTGAATCTCAAAGCTGGGATTATACTACCTATCTTTATACCCACAGAGAGGCTGGCACATGGCAGCCTGTGACAGGAATTCAGAGTAAACTAACATCCGCCATCACACCCGCACTCACTCCCTTCCCTGCCTTCCTCCCACAGGGTGTTTGGCTTCCCCGTCCACTACACAGACGTCTCCAACATGAGCCGCTTGGCGAGGCAGAGACTGCTGGGCCGATCGTGGAGCGTGCCGGTCATCCGCCACCTCTTCGCTCCGCTGAAGGAATATTTTGCTTGTGTGTAAGGGACATGGGGGCAAACTGAAGTAGTGATGATAAAAAAGTTAAACAAACAAACAAACAAAAAACAAAACAAAACAATAAAACACCAAGAACGAGAGGACGGAGAAAAGTTCAGCACCCAGAAGAGAAAAAGGAATTTAAAGCAAACCACAGAGGAGGAAAACGCCGGAGGGCTTGGCCTTGCAAAAGGGTTGGACATCATCTCCTGAGTTTTCAATGTTAACCTTCAGTCCTATCTAAAAAGCAAAATAGGCCCCTCCCCTTCTTCCCCTCCGGTCCTAGGAGGCGAACTTTTTGTTTTCTACTCTTTTTCAGAGGGGTTTTCTGTTTGTTTGGGTTTTTGTTTCTTGCTGTGACTGAAACAAGAGAGTTATTGCAGCAAAATCAGTAACAACAAAAAGTAGAAATGCCTTGGAGAGGAAAGGGAGAGAGGGAAAATTCTATAAAAACTTAAAATATTGGTTTTTTTTTTTTTCCTTTTCTATATATCTCTTTGGTTGTC", [400,600], [[250,500]]))
"""test = generate_grna("ATGCGTGCCTGACCTGTATTTAGTTTTTAGCAGAGCCAGGGCTGCCTTTGGCCCCAGCTGCTGGGGCAGTTAGGCATCTTAGTTATGCATTAGTTTGTTAATTGCTAGCTTGGAAGTTTGGGGTTCCTCAGCTCTGTTCTTCCTACAGGAGTGTGTGAGGGAGAAATTAAACTTACTGCATGAATTCCTGCAAACAGAAATAAAAAGCCAGTTGTGTGACTTGGAAACCAAATTACATAAAGAGGAATTATCTGAGGTAAGTCTGTCCTTTTTCCCAGTTTCCAGAAAGCTACCTGCTTGTGACCAGAGGCAGAGGCCCCTGTGTGTCCTACTGTGATGTCTTCCTACTAAGATATCCTTTCCTGTTGTAGGAAGGCTACCTGGCTAAAGTCAAGTCCCTCTTAAATAAGGATTTGTCCTTGGAGAACGGAACACACACTCTCACTCAAAAAGCCAACGGTTGTCCCGCCAACGGGAGCCGGCCAACCTGGAGAGCAGAAATGGCAGACTCAAATAGATCCCCAAGATCCAGGCCCAAGCCTCGGGGACCCAGGAGAAGCAAGTCGGACAGTGACACCCTTTGTAAGGACACGAGACATACTGGTTTGTTGGTTGTTTTGGGTTTTCTGGTTGTTTTTAAGACGGTACCTGTCAGTATAGGTCAGCCTGGCCTTAGAAGTCTGTCCTGCTTCAGCCCCCTGAGTGCTGGTGTGATAGGCCTATGCTGTCTCACAGGTCTCAGTTGACTTGTTTGCATCTTGGAGGCTGTGGTTCATTTGAAGAGTCAGGCCACTTGAACCATGGCACGTACCACATGTATTCTAGTTGTTTGTTTCTGGTGCTGAGAAGGAAGCCTAAGGCCTCGAACACTAAGCCCATGCTGACACTCAGCTGCGTGTCTGTGCCCCTGCCTCTGGAGTGCTAGTGGTATAACCCAGGACCTGGGTGTGCTAGGGAAGGGCACATCTATCTGCCAAGCCCTGGCTGGCTTTTTCAGGGGAGG", "seqnterm")

for grna in test.keys():
    print(test[grna])"""