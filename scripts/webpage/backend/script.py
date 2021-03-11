# -*- coding: utf-8 -*-
# @Time  : 2021/3/9 1:36 PM
# @Author :
# Apply Python to run R script
import subprocess


def execute_cmd(r_path, r_program, args):
    """
    run R through cmd
    :param r_path: R absolute path
    :param r_program: R programming pah
    :param args: R file path
    :return: {"code": 0, "message": ""}
    """
    cmd = '"{}" "{}" "{}"'.format(r_path, r_program, args)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    res = {'status': True, 'message': 'successful'}
    with open('./output/log.txt', 'w') as f:
        for line in p.stdout.readlines():
            line = line.decode('gbk')        # windows:gbk ios:utf-8
            f.write(line)
            if 'Error:' in line:
                res['status'] = False
                res['message'] = line
    return res



