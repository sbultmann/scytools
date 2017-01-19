from django import forms


class RequestForm(forms.Form):
    CHOICES = (('Mouse', 'Mouse'),('Human', 'Human'),)
    gene_name = forms.CharField(label='gene name', max_length=100)
    organism = forms.CharField(label='organism', widget=forms.Select(choices=CHOICES))

