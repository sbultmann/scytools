from django.db import models

class test(models.Model):
    gene_name = models.CharField(max_length=200)
    org_name = models.CharField(max_length=200)
    gene_org = models.CharField(max_length=200)


