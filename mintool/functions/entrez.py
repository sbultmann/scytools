from Bio import Entrez


Entrez.email = "bultmann@bio.lmu.de"
Entrez.tool = "mintool"

handle = Entrez.esearch(db="gene", term="dnmt3b[Gene Name]) AND Mouse[Organism] ")
record = Entrez.read(handle)
ids = record["IdList"]


handle = Entrez.efetch(db="gene", id=ids, rettype="gb", retmode="xml")
records = Entrez.parse(handle)
#description: ['Entrezgene_summary']
#chromosome: ['Entrezgene_locus'][0]['Gene-commentary_label']

for record in records:
    #print(record['Entrezgene_summary'])
    x = record['Entrezgene_locus'][0]['Gene-commentary_products']
    #print(record['Entrezgene_locus'][0].keys())
    print(record['Entrezgene_locus'][0]['Gene-commentary_seqs'])
    locus_accession = record['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']\
        ['Seq-interval']['Seq-interval_id']['Seq-id']['Seq-id_gi']
    locus_from =record['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']\
        ['Seq-interval']['Seq-interval_to']
    for i in x:
        if 'Gene-commentary_products' in i.keys() and "X" not in i['Gene-commentary_label']:
     #       print(i['Gene-commentary_seqs'])
            isoform_name = i['Gene-commentary_label']
            isoform_accessionP = i['Gene-commentary_products'][0]['Gene-commentary_accession']
            isoform_accessionN = i['Gene-commentary_seqs'][0]['Seq-loc_whole']['Seq-id']['Seq-id_gi']
            isoform_locs = i['Gene-commentary_products'][0]['Gene-commentary_genomic-coords']\
                [0]['Seq-loc_mix']['Seq-loc-mix']

