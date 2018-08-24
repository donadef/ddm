# -*- coding: utf-8 -*-
import os
import shutil
import glob
import numpy as np
import configparser


def check_step(file):
    if not os.path.isfile(file):
        print("""
        ERROR : The process ended cause the file """ + file + """ has not been created as it should.
                Please check the output above to find the error.""")
        exit()


def clean_md_files():
    for f in glob.glob("*.gro"):
        os.remove(f)
    for f in glob.glob("*.tpr"):
        os.remove(f)
    for f in glob.glob("*.trr"):
        os.remove(f)
    for f in glob.glob("*.edr"):
        os.remove(f)
    for f in glob.glob("*.log"):
        os.remove(f)
    for f in glob.glob("*.xtc"):
        os.remove(f)
    for f in glob.glob("*.mdp"):
        os.remove(f)
    for f in glob.glob("*.cpt"):
        os.remove(f)
    for f in glob.glob("*#"):
        os.remove(f)


def clean_tmp():
    if os.path.isfile('TMP'):
        os.remove('TMP')


def compute_std(col, file):
    lcol = []
    with open(file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                c = line.lstrip(' ').rstrip('\n').split(' ')[col - 1]
                lcol.append(float(c))
    a = np.array(lcol)
    # Compute the standard deviation
    std = np.std(a)
    return float(std)


def compute_kf(std):
    return (8.314 * 300 / 1000) / std ** 2


def compute_kf_plus(kf):
    if kf > 10:
        n = np.format_float_scientific(kf)
        kf_plus = str(n.split('e')[0][0])
        for e in range(int(n.split('e')[1]) + 1):
            kf_plus += '0'
    else:
        kf_plus = 100
    return float(kf_plus)

def compute_fluct(x0, lam, kf, col, rms_file, dihe=False):
    sum = 0
    cc = 0
    with open(rms_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                c = line.lstrip(' ').rstrip('\n').split(' ')[col - 1]
                delta = float(c) - x0
                if dihe:
                    absd = abs(delta)
                    if absd > np.pi:
                        delta = (2 * np.pi) - absd
                    else:
                        delta = absd
                sum += delta ** 2
                cc += 1
    return (float(kf) / 2) * (sum / cc)


def compute_trapez(rms_file, verbose=False):
    """ Compute the trapezoidal integration"""
    l = []
    with open(rms_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                ll, dhdl = line.lstrip(' ').rstrip('\n').split(' ')
                l.append([float(ll), float(dhdl)])
    sum = 0
    for i in range(len(l) - 1):
        dA = 0.5 * (l[i][1] + l[i+1][1]) * (l[i][0] + l[i+1][0])
        sum += dA
    return sum / 4.184

class DDMClass:
    def __init__(self, config, complex):
        self.config = config

        self.host = self.config['main']['host']
        self.dest = self.config['main']['dest']

        self.complex = complex
        self.ipdb = os.path.basename(self.complex).rstrip('.pdb')

        self.static_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'static')
        self.awk_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'awk')

        self.directory = self.dest

        self.files_to_store = []

    def store_files(self):
        store_dir = os.path.join(self.directory, 'STORE')
        if not os.path.exists(store_dir):
            os.makedirs(store_dir)

        for file in self.files_to_store:
            try:
                shutil.copy(file, store_dir)
            except FileNotFoundError:
                print('File ' + file + ' has not been created.')
                exit()


class Guest:
    def __init__(self, guest_name, complex, dest):
        self.name = guest_name
        self.complex = complex
        self.dest = dest

        self.pdb_file_path = os.path.join(self.dest, self.name + '.pdb')
        self.pdb = self.create_pdb()

        self.beg = self.find_beg()
        self.end = self.find_end()

    def create_pdb(self):
        pdb = []
        if not os.path.isfile(self.pdb_file_path):
            with open(self.complex, 'r') as file:
                for line in file:
                    if self.name in line:
                        pdb.append(line)
            f = open(self.pdb_file_path, 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', self.pdb)))
            f.close()
        else:
            with open(self.pdb_file_path, 'r') as file:
                for line in file:
                    pdb.append(line.rstrip('\n'))
        return pdb

    def find_beg(self):
        for i in range(len(self.pdb)):
            if 'ATOM' in self.pdb[i]:
                return self.pdb[i].split(' ')[1]

    def find_end(self):
        for i in range(1, len(self.pdb)+1):
            if 'ATOM' in self.pdb[i]:
                return self.pdb[i].split(' ')[1]
