from os import path
from shutil import copyfile
import os

key = " call "
src_dir = '/home/ycmstker/src'
new_dir = path.join(src_dir, 'form')

def run_error(filename, line):
    return "filename: {}, line: {}".format(filename, line)


def dfs(filename):
    fpath = path.join(src_dir, filename)
    assert path.exists(fpath), "filename: {}".format(filename)
    dpath = path.join(new_dir, 'form_'+filename)
    copyfile(fpath, dpath)
    with open(filename, "r") as f:
        lines = f.readlines()
    for line in lines:
        if line[0] in ['*', '!', 'c', 'C']:
            continue
        if key in line:
            cnt = line.count(key)
            assert cnt == 1, run_error(filename, line)
            ind = line.find(key)
            ind2 = ind + 6
            assert line[ind2] != ' ', run_error(filename, line)
            ind3 = line[ind2:].find('(')
            newfilename = line[ind2:ind2+ind3]+'.f'
            if newfilename in ['exit.f']:
                continue
            dfs(newfilename)

if __name__ == '__main__':
    filename = 'e_c3d.f'
    dfs(filename)
