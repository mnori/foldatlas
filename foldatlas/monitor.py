import os

# must call "sudo apt-get install sendmail" first...


# if sts != 0:
#     print("Sendmail exit status "+str(sts))



def send_error(recipient, error_details):
	SENDMAIL = "/usr/sbin/sendmail" # sendmail location

	p = os.popen("%s -t" % SENDMAIL, "w")
	p.write("To: "+recipient+"\n")
	p.write("Subject: FoldAtlas error\n")
	p.write("\n") # blank line separating headers from body
	p.write("Some text\n")
	p.write("some more text\n")
	sts = p.close()
