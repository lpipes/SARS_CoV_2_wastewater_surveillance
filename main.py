from flask import Flask, render_template, request
import atexit, yaml

app = Flask(__name__, template_folder='../frontend/templates', static_folder='../frontend/static')
account_dic = {}


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
            save()
            msg = 'Account created - Login now!'
            return render_template('index.html', msg=msg)

    msg = 'Account created - Login now!'
    return render_template('index.html', msg=msg)


# Method to add new accounts into dictionary
def create_account(username, password):
    account_dic[username] = password


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
def save():
    with open('storage.yml', 'w') as f:
        yaml.dump(account_dic, f)


# Load yml file into dictionary
def load():

    with open('storage.yml', 'r') as file_descriptor:
        data = yaml.safe_load(file_descriptor)

        for key in data:
            account_dic[key] = data.get(key)


if __name__ == "__main__":
    load()
    app.run()
