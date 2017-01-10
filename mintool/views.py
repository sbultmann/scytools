from django.shortcuts import render

from django.http import HttpResponse, HttpResponseRedirect
from .forms import RequestForm
from .models import test
from .functions.entrez import *
import json


def index(request):
    if request.method == 'POST':
        # create a form instance and populate it with data from the request:
        form = RequestForm(request.POST)
        # check whether it's valid:
        if form.is_valid():
            gene_name = form.cleaned_data['gene_name']
            organism = form.cleaned_data['organism']
            gene_org = gene_name + organism
            request = test.objects.create(
                gene_name=gene_name,
                org_name=organism,
                gene_org=gene_org, )
            return HttpResponseRedirect('%s/request/' % str(request.id))
    else:
        form = RequestForm()
    return render(request, 'mintool/index.html', {'form': form})


def requestNCBI(request, gene_org):
    info = test.objects.get(pk=gene_org)
    res = query_ncbi(str(info.gene_name), str(info.org_name))
    return render(request, 'mintool/main_results.html', {'result': res})


def gRNAdesign(request):
    if request.method == 'GET':
        from Bio.Seq import Seq
        sequence = request.GET['sequence']
        term = request.GET['term']

        response_data = dict()
        response_data["gRNAs"] = generate_grna(sequence, term)
        response_data["primers"] = design_primers(sequence, [200, 500], [[250, 450]])
        response_data["sequence"] = sequence
        response_data["rev_sequence"] = str(Seq(sequence).complement())
        response_data["term"] = term

        return HttpResponse(
            json.dumps(response_data),
            content_type="application/json"
        )
    else:
        return HttpResponse(
            json.dumps({"nothing to see": "this isn't happening"}),
            content_type="application/json"
        )


def create_csv(request):
    if request.method == 'GET':
        toligo = request.GET['toligo']
        gRNA_seq = request.GET['gRNA-seq']
        primer_f = request.GET['primer_f']
        primer_r = request.GET['primer_r']
        f_oligo, r_oligo = gRNA_oligos(gRNA_seq)
        locus_sequence = request.GET['locus_sequence']
        csv = "targeting_oligo,"+toligo+"\n"
        csv = csv + "screening_primer_f," + primer_f + "\n"
        csv = csv + "screening_primer_r," + primer_r + "\n"
        csv = csv + "gRNA_upper," + f_oligo + "\n"
        csv = csv + "gRNA_lower," + r_oligo + "\n"
        csv = csv + "full_locus_sequence," + locus_sequence + "\n"
        response = HttpResponse(csv, content_type='application/pdf')
        response['Content-Disposition'] = 'attachment; filename="targeing_oligos.csv"'
        return response
    else:
        csv = "Nothing to see - something went wrong"
        response = HttpResponse(csv, content_type='application/pdf')
        response['Content-Disposition'] = 'attachment; filename="targeing_oligos.csv"'
        return response
    #toligo, gRNA-seq, primer_f, primer_r

