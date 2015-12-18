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
			print("It worked!")
		else:
			send_error(text)
	except:
		send_error(traceback.format_exc())

def send_error(error_details):
	print("FAILED")
	SENDMAIL = "/usr/sbin/sendmail" # sendmail location

	p = os.popen("%s -t" % SENDMAIL, "w")
	p.write("To: "+recipient+"\n")
	p.write("Subject: FoldAtlas error\n")
	p.write("\n") # blank line separating headers from body
	p.write(error_details)
	sts = p.close()

run_test()