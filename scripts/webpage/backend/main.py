import re

from flask import Flask, render_template, request, flash, redirect, url_for, make_response
import atexit, yaml
import os
from werkzeug.utils import secure_filename
from flask import send_from_directory
from flask_bootstrap import Bootstrap
from flask_mail import Mail,Message  # sending the emails

from pic_plot.plot_pic import picture
from script import execute_cmd
from send_email import send_email

basedir = os.path.abspath(os.path.dirname(__file__))

uploadDir = os.path.join(basedir, 'upload')

# Rpath, with packages installed
'''
Note: Please enter your Rscript path here!!!
'''
Rscript_PATH = 'D:/R/R-4.0.4/bin/Rscript'
# output path
OUTPUT_PATH = os.path.join(basedir, 'output')

# set the path for uploaded files
UPLOAD_FOLDER = uploadDir

# set the form of files
ALLOWED_EXTENSIONS = {'txt'}

app = Flask(__name__, template_folder='../frontend/templates', static_folder='../frontend/static')

# email address information
app.config.update(dict(
    DEBUG=True,
    MAIL_SERVER="smtp.163.com",
    # MAIL_PORT=587,
    # MAIL_USE_TLS=True,
    MAIL_PORT=465,  # 163--465
    MAIL_USE_SSL=True,
    MAIL_USE_TLS=False,
    MAIL_USERNAME="sars_cov2@163.com",   # my email address
    MAIL_PASSWORD="BZIXVTMTZUDTGVWL",          # my email password
    MAIL_DEFAULT_SENDER=("sars_cov2@163.com"),  # defult sending address
    MAIL_DEBUG=True,

))


mail = Mail(app)   # create a mail

bootstrap = Bootstrap(app)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['SECRET_KEY'] = os.urandom(24)
account_dic = {}
email_list = []
email_file = {}


@app.route("/")
def main():
    return render_template("index.html")

@app.route("/login", methods=['GET', 'POST'])
def login():

    if request.method == 'POST':

        # Requesting from web page
        username = request.form['username']
        password = request.form['password']

        # Making sure account exists
        if has_account(username.lower()):

            # Making sure password matches
            if password_validator(username.lower(), password):
                return render_template("end_point.html")

            else:
                error = 'Incorrect password'
                return render_template('index.html', error=error)

        else:
            error = 'Username is not found'
            return render_template('index.html', error=error)

    return render_template("end_point.html")


@app.route("/register")
def register():
    return render_template("register.html")


@app.route("/register-handler", methods=['GET', 'POST'])
def register_handler():

    if request.method == 'POST':

        # Requesting from web page
        username = request.form['username']
        password = request.form['password']
        re_password = request.form['re-password']

        # A bunch of small checks
        if 4 > len(password) < 20:
            error = 'Passwords must be 4 - 20 characters.'
            return render_template('register.html', error=error)

        elif 4 > len(username) < 20:
            error = 'Passwords must be 4 - 20 characters.'
            return render_template('register.html', error=error)

        elif password != re_password:
            error = 'Passwords do not match.'
            return render_template('register.html', error=error)

        elif has_account(username.lower()):
            error = 'Username is not available.'
            return render_template('register.html', error=error)

        else:
            create_account(username.lower(), password)
            save_account()
            msg = 'Account created - Login now!'
            return render_template('index.html', msg=msg)

    msg = 'Account created - Login now!'
    return render_template('index.html', msg=msg)


@app.route("/get_email", methods=['GET', 'POST'])
def get_email():
    if request.method == 'POST':
        email = request.form['email']
        if email == '':
            error = 'Please enter your email address!'
            return redirect('end_point.html',error=error)
        else:
            email_list.append(email)
            print('email_list:',email_list)
            return render_template('upload.html')


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/upload', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        # check if the post request has the file part
        print(request.files)
        if 'file' not in request.files:
            # flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        # check if the path exist, if not, create it
        if not os.path.exists(uploadDir):
            os.makedirs(uploadDir)
        if file:
            if allowed_file(file.filename):
                filename = secure_filename(file.filename)
                upload_file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                file.save(upload_file_path)
                print('----:',email_list, filename)

                # import R, process the data and send emails
                res = execute_cmd(Rscript_PATH, 'EM.R', upload_file_path)
                if not res['status']:
                    error_msg = re.search('Error: (.*?)"', res['message']).group(1)
                    send_email(filelist=[], pic_path='', receiver=email_list[-1], content=error_msg)
                    return render_template('end_point.html', msg=error_msg)
                create_email(email_list, filename)
                save_email()
                 # Use R for csv and picture results
                csv_path = os.path.join(OUTPUT_PATH, 'r_output.csv')
                pdf_path = os.path.join(OUTPUT_PATH, 'r_output.pdf')
                picture(csv_path=csv_path, out_path=OUTPUT_PATH)
                pic_path = os.path.join(OUTPUT_PATH, 'pic_plot.png')
                send_email(filelist=[pdf_path, csv_path], pic_path=pic_path, receiver=email_list[-1])

                msg = 'Upload Load Successful!'
                print(msg)
                return render_template('end_point.html', msg=msg)
            else:
                flash('Unknown Types!', 'danger')
                return render_template('end_point.html')
        else:
            flash('No File Selected.', 'danger')
        return render_template('end_point.html')
    return render_template('end_point.html')

@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename)

# Method to add new accounts into dictionary
def create_account(username, password):
    account_dic[username] = password

# Mewhod to add new email addresses into dictionary
def create_email(email_list, filename):
    email_file[email_list[-1]] = filename

# Method to check if password is valid
def password_validator(username, password):
    if account_dic is None:
        return False
    check_password = account_dic[username]
    if check_password == password:
        return True
    return False

# Method to check for existing account
def has_account(username):
    if account_dic is None:
        return False
    for user in account_dic.keys():
        if user == username:
            return True
    return False

# Save dictionary into a yml
def save_account():
    with open('storage.yml', 'w') as f:
        yaml.dump(account_dic, f)

# Save dictionary into a yml
def save_email():
    with open('email_file.yml', 'a') as f2:
    	yaml.dump(email_file,f2)

# Load yml file into dictionary
def load():
    with open('storage.yml', 'r') as file_descriptor:
        data = yaml.safe_load(file_descriptor)
        for key in data:
            account_dic[key] = data.get(key)

# from main import send_email_by_auto
def send_email_by_auto():
    import requests
    respons = requests.get('127.0.0.1:5000/email_send_attach')
    print(respons.status_code)  #return 200 is success


if __name__ == "__main__":
    load()
    app.run()
