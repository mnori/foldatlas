from django.conf.urls import patterns, include, url
from django.contrib import admin
import rnabrowserapp.views as views

urlpatterns = patterns('',
	url(r'^admin/', admin.site.urls),
	url(r'^$', views.index),
	url(r'^transcript/([a-zA-Z0-9]+)/$', views.transcript)
)
