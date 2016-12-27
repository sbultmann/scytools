from django import forms


class RequestForm(forms.Form):
    CHOICES = (('mm', 'Mouse'),('hs', 'Human'),)
    gene_name = forms.CharField(label='gene name', max_length=100)
    organism = forms.CharField(label='organism',widget=forms.Select(choices=CHOICES))
