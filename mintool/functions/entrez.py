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
                    loc_to = int(loc['Seq-loc_int']['Seq-interval']['Seq-interval_to'])
                    loc_from = int(loc['Seq-loc_int']['Seq-interval']['Seq-interval_from'])
                    real_exons.append(loc_to)
                    real_exons.append(loc_from)
                    loc_to = abs(int((float(loc['Seq-loc_int']['Seq-interval']['Seq-interval_to'])-maxlenght)*ratio))
                    loc_from = abs(int((float(loc['Seq-loc_int']['Seq-interval']['Seq-interval_from'])-maxlenght)*ratio))

                    loc_strand = str(loc['Seq-loc_int']['Seq-interval']['Seq-interval_strand']).split(":")[2].split("}")[0][2:-1]
                    if loc_to > loc_from:
                        exon = [loc_to, loc_from]
                    else:
                        exon = [loc_from, loc_to]
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


#print(query_ncbi("Dnmt3a", "Mouse"))

#print(extract_sequences(372099101,20907520,20941453,2))
#print(extract_sequences(372099098,3835197,3907746,1))

#20908751 20908471 2 372099101