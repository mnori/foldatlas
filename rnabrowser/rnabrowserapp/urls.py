from django.conf.urls import patterns, url

from rnabrowserapp import views

urlpatterns = patterns('',
    url(r'^$', views.index, name='index'),
)
