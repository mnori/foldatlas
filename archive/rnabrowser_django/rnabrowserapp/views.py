# what does this do?
# from django.shortcuts import render

from django.http import HttpResponse
from django.template import RequestContext, loader
from rnabrowserapp.models import Transcript, TranscriptSequence
from django.db.models import Q
from django.shortcuts import redirect

# Define index view.
def index(request):

    error = ""
    query = request.POST.get("query")
    if query:
        print("Query: "+query)

        try:
            # if transcript exists, redirect to its page.
            transcript = Transcript.objects.get(Q(id=query)|Q(gi_number=query))

            # print(transcript)

            # Redirects to whatever Transcript.get_absolute_url() is
            return redirect(transcript)

            # do something with transcript.id

        except Transcript.DoesNotExist:
            # if it doesn't exist, this will show an error underneath the search box.
            error = "No sequence found matching \""+str(query)+"\""

    template = loader.get_template('rnabrowserapp/index.html')
    context = RequestContext(request, {
        'error': error,
    })
    return HttpResponse(template.render(context))

def transcript(request, transcript_id):

    transcript = Transcript.objects.get(id=transcript_id)
    transcript_sequence = TranscriptSequence.objects.get(transcript_id=transcript_id)

    # get the transcript using its ID

    # if transcript was not found, redirect to the home page

    # otherwise display the sequence

    template = loader.get_template('rnabrowserapp/transcript.html')
    context = RequestContext(request, {
        'transcript': transcript,
        'transcript_sequence': transcript_sequence
    })
    return HttpResponse(template.render(context))