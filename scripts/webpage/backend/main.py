from flask import Flask, render_template, request, flash, redirect, url_for
import atexit, yaml
import os
from werkzeug.utils import secure_filename
from flask import send_from_directory

basedir = os.path.abspath(os.path.dirname(__file__))
uploadDir = os.path.join(basedir, 'upload')

# set the path for uploaded files
UPLOAD_FOLDER = uploadDir

# set the form of files
ALLOWED_EXTENSIONS = {'txt', 'fasta'}


app = Flask(__name__, template_folder='../frontend/templates', static_folder='../frontend/static')
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
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
            return render_template('upload.html')

    
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        # check if the post request has the file part
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        # if user does not select file, browser also
        # submit an empty part without filename
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            save_email(email_list,filename)
            return redirect(url_for('uploaded_file',
                                    filename=filename))
    return '''
    <!doctype html>
    <title>Upload new File</title>
    <h1>Upload new File</h1>
    <form method=post enctype=multipart/form-data>
      <input type=file name=file>
      <input type=submit value=Upload>
    </form>
    '''

@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'],
                               filename)

# Method to add new accounts into dictionary
def create_account(username, password):
    account_dic[username] = password

#Method to add new emails and fastas into dictionary
# def create_address(email, filename):
# 	email_file[email] = filename

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
def save_email(email_list,filename):
    for email in email_list:
        
        email_file[email] = filename

    with open('email_file.yml', 'w') as f2:

    	yaml.dump(email_file,f2)

# Load yml file into dictionary
def load():

    with open('storage.yml', 'r') as file_descriptor:
        data = yaml.safe_load(file_descriptor)

        for key in data:
            account_dic[key] = data.get(key)

if __name__ == "__main__":
    load()
    app.run()
