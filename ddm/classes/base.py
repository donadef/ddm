# -*- coding: utf-8 -*-
import os
import shutil
import glob
import numpy as np
from scipy import constants as c
import math


ORGANIZE = {'modeling': 'modeling',
            'solvate-bound': 'solvate-bound',
            'solvate-unbound': 'solvate-unbound',
            'pick-reference': 'pick-reference',
            'monitor-CVs': 'monitor-CVs',
            'confine-bound': 'confine-bound',
            'vba-bound': 'vba-bound',
            'alchemical-bound': 'alchemical-bound',
            'confine-unbound': 'confine-unbound',
            'vba-unbound': 'vba-unbound',
            'alchemical-unbound': 'alchemical-unbound',
            }


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


def compute_mean(col, file):
    lcol = []
    with open(file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                c = line.lstrip(' ').rstrip('\n').split(' ')[col - 1]
                lcol.append(float(c))
    a = np.array(lcol)
    # Compute the standard deviation
    mean = np.mean(a)
    return float(mean)


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


def compute_kf(std, temp):
    return (8.314 * temp / 1000) / std ** 2


def compute_kf_plus(kf):
    if kf > 10:
        n = np.format_float_scientific(kf)
        kf_plus = str(n.split('e')[0][0])
        for e in range(int(n.split('e')[1]) + 1):
            kf_plus += '0'
    else:
        kf_plus = 100
    return float(kf_plus)


def compute_fluct(x0, kf, col, file, dihe=False):
    sum = 0
    cc = 0
    file_content = open(file, 'r')
    file_lines = file_content.readlines()
    file_content.close()
    perc = 1
    list_length = int(len(file_lines) * perc)
    for i in range(list_length):
        line = file_lines[i]
        if not line.startswith('#'):
            c = line.lstrip(' ').rstrip('\n').split(' ')[col - 1]
            delta = float(c) - float(x0)
            if dihe:
                absd = abs(delta)
                if absd > np.pi:
                    delta = (2 * np.pi) - absd
                else:
                    delta = absd
            sum += delta ** 2
            cc += 1
    return (float(kf) / 2) * (sum / cc)


def compute_trapez(fluct_list, col):
    """ Compute the trapezoidal integration"""
    l = []
    for item in fluct_list:
        item = item.lstrip(' ').rstrip('\n').split(' ')
        ll = item[0]
        dhdl = item[col - 1]
        l.append([float(ll), float(dhdl)])
    sum = 0
    for i in range(len(l) - 1):
        dA = 0.5 * (l[i][1] + l[i+1][1]) * (l[i+1][0] - l[i][0])
        sum += dA
    return sum / 4.184


def compute_work(kappa, temp):
    kT = c.Boltzmann * c.N_A * temp / 1000  # kJ/mol

    #  positional restrain
    Ztr = 1.66058  # nm^3/molecule (volume/molecule @ 1M)
    Ztr_R = (kappa[6] ** 2) * math.sin(kappa[7]) * ((2 * np.pi * kT) ** 1.5) / math.sqrt(kappa[0] * kappa[1] * kappa[2])

    # orientation restrain
    Zrot = 8 * (np.pi ** 2)
    Zrot_R = math.sin(kappa[8]) * ((2 * np.pi * kT) ** 1.5) / math.sqrt(kappa[3] * kappa[4] * kappa[5])

    # Restraint work
    Wr = -kT * (math.log(Ztr/Ztr_R) + math.log(Zrot/Zrot_R)) / 4.184  # kcal/mol
    return Wr


class DDMClass:
    def __init__(self, config):
        self.config = config

        self.dest = self.config['main']['dest']
        self.ff_param = self.config['main'].get('ff_parameters', False)
        self.temp = float(self.config['main'].get('temperature', '300'))

        self.static_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'static')
        self.awk_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'awk')

        self.directory = self.dest

    def run(self):
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

        os.chdir(self.directory)

    def store_files(self, files):
        store_dir = os.path.join(self.directory, 'STORE')
        if not os.path.exists(store_dir):
            os.makedirs(store_dir)

        for file in files:
            try:
                shutil.copy(file, store_dir)
            except FileNotFoundError:
                print('File ' + file + ' has not been created.')
                exit()


class Guest:
    def __init__(self, guest_name, complex_obj, dest):
        self.name = guest_name
        self.complex = complex_obj
        self.dest = dest

        self.pdb_file_path = os.path.join(self.dest, self.name + '.pdb')
        self.pdb = self.create_pdb()

        self.beg = self.find_beg()
        self.nb_at = self.find_nb_atom()
        self.end = self.find_end()

    def create_pdb(self):
        pdb = []
        if not os.path.isfile(self.pdb_file_path):
            for line in self.complex.pdb:
                if self.name in line:
                    pdb.append(line)
            f = open(self.pdb_file_path, 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', pdb)))
            f.close()
        else:
            with open(self.pdb_file_path, 'r') as file:
                for line in file:
                    pdb.append(line.rstrip('\n'))
        return pdb

    def find_beg(self):
        for i in range(len(self.pdb)):
            if 'ATOM' in self.pdb[i]:
                return self.pdb[i].split()[1]

    def find_end(self):
        for i in range(1, len(self.pdb)+1):
            if 'ATOM' in self.pdb[-i]:
                return self.pdb[-i].split()[1]

    def find_nb_atom(self):
        nb = 0
        for i in range(len(self.pdb)):
            if 'ATOM' in self.pdb[i]:
                nb += 1
        return nb

class Host:
    def __init__(self, host_name, complex_obj, dest):
        self.name = host_name
        self.complex = complex_obj
        self.dest = dest

        self.pdb_file_path = os.path.join(self.dest, self.name + '.pdb')
        self.pdb = self.create_pdb()

    def create_pdb(self):
        pdb = []
        if not os.path.isfile(self.pdb_file_path):
            for line in self.complex.pdb:
                if self.name in line:
                    pdb.append(line)
            f = open(self.pdb_file_path, 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', pdb)))
            f.close()
        else:
            with open(self.pdb_file_path, 'r') as file:
                for line in file:
                    pdb.append(line.rstrip('\n'))
        return pdb


class Complex:
    def __init__(self, pdb_complex, dest):
        self.pdb_ori = pdb_complex
        self.name = os.path.basename(self.pdb_ori).rstrip('.pdb')
        self.dest = dest

        self.pdb_file_path = os.path.join(self.dest, os.path.basename(self.pdb_ori))

        self.pdb = self.create_pdb()

    def create_pdb(self):
        pdb = []
        if not os.path.isfile(self.pdb_file_path):
            shutil.copy(self.pdb_ori, self.dest)
        with open(self.pdb_file_path, 'r') as file:
            for line in file:
                pdb.append(line.rstrip('\n'))
        return pdb
