from django.shortcuts import render

from django.http import HttpResponse,HttpResponseRedirect
from .forms import RequestForm
from .models import test


def index(request):
    if request.method == 'POST':
        # create a form instance and populate it with data from the request:
        form = RequestForm(request.POST)
        # check whether it's valid:
        if form.is_valid():
            gene_name = form.cleaned_data['gene_name']
            organism = form.cleaned_data['organism']
            gene_org = gene_name+organism
            request = test.objects.create(
                gene_name=gene_name,
                org_name=organism,
                gene_org=gene_org,)
            return HttpResponseRedirect('%s/request/' % str(request.id))
    else:
        form = RequestForm()
    return render(request, 'mintool/index.html', {'form': form})


def requestNCBI(request, gene_org):
    response = "You're looking at the results of question %s."
    return HttpResponse(response % gene_org)
