from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^mintool/', views.index, name='index'),
    url(r'^(?P<gene_org>[0-9A-Za-z]+)/request/$', views.requestNCBI, name="request"),
    url(r'^[0-9A-Za-z]+/request/gRNA-design/$', views.gRNAdesign, name="gRNAdesign"),
    url(r'^[0-9A-Za-z]+/request/create_csv/$', views.create_csv, name="create_csv"),
]