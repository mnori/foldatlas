# Add this to crontab: 0,30 * * * * python3 /home/ubuntu/foldatlas/foldatlas/monitor.py

import traceback
import os
import urllib.request  # the lib that handles the url stuff

test_url = "http://www.foldatlas.com/transcript/AT2G45180.1"
recipient = "matthew.gs.norris@gmail.com"
search_str = "AT2G45180.1"

def run_test():
    try:
        data = urllib.request.urlopen(test_url) # it's a file like object and works just like a file
        text = str(data.read())

        if search_str in text:
            send_mail("FoldAtlas success", "It worked!")
            print("It worked!")
        else:
            send_mail("FoldAtlas error", text)
    except:
            send_mail("FoldAtlas error", traceback.format_exc())

def send_mail(subject, body):
    SENDMAIL = "/usr/sbin/sendmail" # sendmail location
    
    p = os.popen("%s -t" % SENDMAIL, "w")
    p.write("To: "+recipient+"\n")
    p.write("From: FoldAtlas\n")
    p.write("Subject: "+subject+"\n")
    p.write("\n") # blank line separating headers from body
    p.write(body)
    sts = p.close()

run_test()
