import smtplib
import email.mime.text
import email.mime.multipart
import yaml

import os
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.application import MIMEApplication
from email.mime.image import MIMEImage

basedir = os.path.abspath(os.path.dirname(__file__))
path = os.path.join(basedir, 'upload')

key = '.txt'
def get_type_file(path, key): # Get the output txt files
	# list all the filenames
    files = os.listdir(path)
    file_list=[]
    for file in files:
        if key in file:
            file_list.append(file)
    return file_list

def send_email(filelist, pic_path, receiver, content=None):
	load()
	smtpHost = 'smtp.163.com'
	sendAddr = 'sars_cov2@163.com'
	password = 'YTKPXXRKISOTASSV'
	subject = "SARS_COV2_WasteWater_Surveillance_Result"
	content = content if content else "Here's the result."
	receiver = receiver #lpipes@berkeley.edu

	msg = MIMEMultipart()

	msg['from'] = sendAddr
	msg['to'] = receiver
	msg['Subject'] = subject
	txt = MIMEText(content, 'plain', 'utf-8')
	msg.attach(txt)  # Input the main body

	for file in filelist:  # Add the attachment(output files)
		filename = file
		part = MIMEApplication(open(filename, 'rb').read())
		part.add_header('Content-Disposition', 'attachment', filename=filename)
		msg.attach(part)
	if pic_path:
		figure = MIMEApplication(open(pic_path, 'rb').read())
		figure.add_header('Content-Disposition', 'attachment', filename='outcome.png')
		msg.attach(figure)

	server = smtplib.SMTP(smtpHost, 25)
    # server.set_debuglevel(1)  # can be used to check bugs
	server.login(sendAddr, password)
	server.sendmail(sendAddr, receiver, str(msg))
	   # print("\n"+ str(len(filelist)) + "files are successfully sent!")
	server.quit()

email_file = {}
# Load email addresses and files
def load():
    with open('email_file.yml', 'r') as file_descriptor:
        data = yaml.safe_load(file_descriptor)
        for key in data:
            email_file[key] = data.get(key)

def main():
    file_list = get_type_file(path, key)
    send_email(file_list)
 

if __name__ == '__main__':
	main()
