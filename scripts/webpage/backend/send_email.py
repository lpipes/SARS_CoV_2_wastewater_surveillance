import smtplib
import email.mime.text
import email.mime.multipart
import yaml

import os
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.application import MIMEApplication


path = 'C:/Users/61427/Desktop/login/outputs'
key = '.txt'
def get_type_file(path, key): # Get the out put txt files
	# list all the filenames
    files = os.listdir(path)
    file_list=[]
    for file in files:
        if key in file:
            file_list.append(file)
    return file_list

def send_email(filelist):
	load()
	smtpHost = 'smtp.163.com'
	sendAddr = 'sars_cov2@163.com'
	password = 'UXRXEQZJIZXLVBPL'
	subject = "SARS_COV2_WasteWater_Surveillance_Result"
	content = "Here's the result."
	receiver = 'guxiaoyi1809@163.com'

	msg = MIMEMultipart()

	msg['from'] = sendAddr
	msg['to'] = receiver
	msg['Subject'] = subject
	txt = MIMEText(content, 'plain', 'utf-8')
	msg.attach(txt)  # Input the main body
	for file in filelist: # Add the attachment(output files)
		filename = file
		part = MIMEApplication(open(path + '/' + filename, 'rb').read())
		part.add_header('Content-Disposition', 'attachment', filename=filename)
		msg.attach(part)

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
 
