# -*- coding: utf-8 -*-
import os
import time
import hashlib
import random
import shutil

from werkzeug.utils import secure_filename
from flask import Flask, render_template, redirect, url_for, request, flash
from flask_uploads import UploadSet, configure_uploads, FASTAS, patch_request_class
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired, FileAllowed
from wtforms import SubmitField

app = Flask(__name__,template_folder='templates', static_folder='static')
# app = Flask(__name__, template_folder='../templates', static_folder='../static')
# set the uploaded file to xiaoyi's path
uploaded_path = '/var/www/FlaskApp/FlaskApp'
app.config['SECRET_KEY'] = 'Waste Water Surveillance'
app.config['UPLOADED_FASTAS_DEST'] = os.path.join(uploaded_path, 'uploaded_file')
# app.config['UPLOADED_FASTAS_DEST'] = os.getcwd() + '/uploaded_file'
app.config['UPLOADED_FASTAS_URL'] = os.path.join(uploaded_path, 'uploaded_file')
# app.config['UPLOADED_FASTAS_URL'] = os.getcwd() + '/uploaded_file'
app.config['UPLOADED_FASTAS_ALLOW'] = set(['fas'])

fastas = UploadSet('fastas', FASTAS)
configure_uploads(app, fastas)
# set the base directory
base_dir = os.path.join(uploaded_path, 'uploaded_file')
# base_dir = os.getcwd() + '/uploaded_file'
user_dir = ''
cmd = ''
paired_reads = ''

class UploadForm(FlaskForm):
    fasta = FileField(validators=[FileAllowed(fastas, u'fastas/txts/fastqs Only!'), FileRequired(u'Choose a file!')])
    submit = SubmitField(u'Upload')


@app.route('/', methods=['GET', 'POST'])
def test():
    if request.method == 'POST':
        return render_template('index.html')
    return render_template('index.html')

@app.route('/uploaded_yes', methods=['POST','GET'])
def upload_yes():
    global user_dir
    global cmd
    paired_reads = 'yes'
    if request.method == 'POST':
        email = request.form['email']
        frequency = request.form['frequency']
        EM_error = request.form['EM_error']
        print(email, frequency, EM_error, paired_reads) # get all the needed info

        # set the number of files
        count = random.randint(100000,1000000)
        # set the user's own file name, e.g. guxiaoyi1809@163.com_1111111
        user_file = email + '_' + str(count)
        # set the file directory for a user
        user_dir = os.path.join(os.path.join(uploaded_path, 'uploaded_file'), user_file)
        # user_dir = os.path.join(os.getcwd() + '/uploaded_file', user_file)
        # if the directory doesn't exist, create one
        if os.path.exists(user_dir):
            user_file = email + '_' + str(count + 1)
            user_dir = os.path.join(os.path.join(uploaded_path, 'uploaded_file'), user_file)
            # user_dir = os.path.join(os.getcwd() + '/uploaded_file', user_file)
            os.mkdir(user_dir)
        else:
            os.mkdir(user_dir)

    form = UploadForm()
    # get the valid files
    if form.validate_on_submit():
        # get all the files the user uploaded
        files_list = request.files.getlist('fasta')
        # print(files_list)
        # for filename in request.files.getlist('fasta'):
        # when the user upload two files and chooses yes
        forward_file = files_list[0]
        reverse_file = files_list[-1]
        # rename the two files the user uploaded
        forward_name = "forward_" + forward_file.filename
        fastas.save(forward_file, name=forward_name)
        reverse_name = "reverse_" + reverse_file.filename
        fastas.save(reverse_file, name=reverse_name)
        # copy the two files to the user's directory
        forward_dir = os.path.join(base_dir, forward_name)
        reverse_dir = os.path.join(base_dir, reverse_name)
        shutil.copy(forward_dir, user_dir)
        shutil.copy(reverse_dir, user_dir)
        # set the forward and reverse path for cmd
        forward_cmd = os.path.join(user_dir, forward_name)
        reverse_cmd = os.path.join(user_dir, reverse_name)

        alignment_path = os.path.join(user_dir, 'alignment.sam')
            # define the path for mismatch.txt
        mismatch_path = os.path.join(user_dir, 'mismatch.txt')
            # turn to the index page before run the command
        cmd = "/home/software/eliminate_strains/./eliminate_strains -i /home/database/msa_0223_updated.fasta.gz " + \
                    "-r /home/database/reference_positions.txt -v /home/database/variants.txt " + \
                    "-e " + EM_error + " " +\
                    "-f " + frequency + " " +\
                    "-o " + mismatch_path + " "\
                    "-s "+ alignment_path + " "\
                    "-d /home/lenore/wuhCor1 -p -1 " + forward_cmd + ' '\
                    "-2 " + reverse_cmd + ' '
        # print(cmd)
        success = True
    else:
        success = False
    return render_template('index_test_yes.html', form=form, success=success)

@app.route('/uploaded_no', methods=['POST','GET'])
def upload_no():
    global cmd
    global user_dir
    paired_reads = 'no'
    if request.method == 'POST':
        email = request.form['email']
        frequency = request.form['frequency']
        EM_error = request.form['EM_error']
        print(email, frequency, EM_error, paired_reads) # get all the needed info

        # set the number of files
        count = random.randint(100000,1000000)
        # set the user's own file name, e.g. guxiaoyi1809@163.com_1111111
        user_file = email + '_' + str(count)
        # set the file directory for a user
        user_dir = os.path.join(os.path.join(uploaded_path, 'uploaded_file'), user_file)
        # user_dir = os.path.join(os.getcwd() + '/uploaded_file', user_file)
        # if the directory doesn't exist, create one
        if os.path.exists(user_dir):
            user_file = email + '_' + str(count + 1)
            user_dir = os.path.join(os.path.join(uploaded_path, 'uploaded_file'), user_file)
            # user_dir = os.path.join(os.getcwd() + '/uploaded_file', user_file)
            os.mkdir(user_dir)
        else:
            os.mkdir(user_dir)

    form = UploadForm()
    # get the valid files
    if form.validate_on_submit():
        files_list = request.files.getlist('fasta')
        filename = files_list[0]
        name = filename.filename
        fastas.save(filename, name=name)
        file_dir = os.path.join(base_dir, name)
        shutil.copy(file_dir, user_dir)
        # define the path for alignment.sam
        alignment_path = os.path.join(user_dir, 'alignment.sam')
        # define the path for mismatch.txt
        mismatch_path = os.path.join(user_dir, 'mismatch.txt')
        # turn to the index page before run the command

        file_cmd = os.path.join(user_dir, name)
        cmd = "/home/software/eliminate_strains/./eliminate_strains -i /home/database/msa_0223_updated.fasta.gz " + \
                    "-r /home/database/reference_positions.txt -v /home/database/variants.txt " + \
                    "-e " + EM_error + " " +\
                    "-f " + frequency + " " +\
                    "-o " + mismatch_path + " "\
                    "-s "+ alignment_path + " "\
                    "-d /home/lenore/wuhCor1 -0 " + file_cmd
        # print(cmd)
        success = True
    else:
        success = False
    return render_template('index_test_no.html', form=form, success=success)

@app.route('/manage')
def manage_file():
    # files_list = os.listdir(app.config['UPLOADED_FASTAS_DEST'])
    files_list = os.listdir(user_dir)
    return render_template('manage.html', files_list=files_list)


@app.route('/open/<filename>')
def open_file(filename):
    file_url = fastas.url(filename)
    
    temp_url = file_url
    path = temp_url[:-len(filename)]

    file_url = os.path.join(path, filename)
    # print(file_url)

    with open(file_url, 'r') as f:
        content = f.read()
    return render_template('browser.html', file_url=file_url, content=content)


@app.route('/delete/<filename>')
def delete_file(filename):
    global user_dir
    file_path = fastas.path(filename)
    os.remove(file_path)
    os.remove(os.path.join(user_dir,filename))
    return redirect(url_for('manage_file'))

@app.route('/submit', methods=['GET', 'POST'])
def submit_file():
    print(cmd)
    return render_template('submit.html', cmd = cmd)

@app.route('/back',methods=['GET', 'POST'])
def back_file():
    os.system(cmd)
    return render_template('submit.html', cmd = cmd)

if __name__ == '__main__':
    app.run('0.0.0.0', port=5001,debug=True)

