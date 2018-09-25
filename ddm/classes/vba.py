# -*- coding: utf-8 -*-
import os
import shutil
import subprocess
import configparser
from scipy import constants as c
import math

from .base import DDMClass, check_step, clean_md_files, compute_mean, compute_fluct, compute_trapez, compute_work, ORGANIZE


class Bound(DDMClass):
    def __init__(self, config, guest, x0, kappa_max, krms_max):
        super(Bound, self).__init__(config)

        self.guest = guest
        self.x0 = x0
        self.kappa_max = kappa_max
        self.krms_max = krms_max

        self.dihe = [False, False, True, False, True, True]

    def create_plumed_tmpl(self):
        f = open(os.path.join(self.static_dir, 'plumed.tmpl'), 'r')
        filedata = f.read()
        f.close()

        newdata = filedata.replace('XXXXX', self.guest.beg)
        newdata = newdata.replace('YYYYY', self.guest.end)
        newdata = newdata.replace('RRRRR', self.guest.name)
        newdata = newdata.replace('REF1', str(self.x0[0]))
        newdata = newdata.replace('REF2', str(self.x0[1]))
        newdata = newdata.replace('REF3', str(self.x0[2]))
        newdata = newdata.replace('REF4', str(self.x0[3]))
        newdata = newdata.replace('REF5', str(self.x0[4]))
        newdata = newdata.replace('REF6', str(self.x0[5]))
        newdata = newdata.replace('KRMS', str(self.krms_max))

        return newdata


class VbaBound(Bound):
    def __init__(self, config, guest, x0, kappa_max, krms_max):
        super(VbaBound, self).__init__(config, guest, x0, kappa_max, krms_max)
        self.directory = os.path.join(self.dest, ORGANIZE['vba-bound'])
        self.static_dir = os.path.join(self.static_dir, 'vba-bound')
        self.prev_store = os.path.join(self.dest, ORGANIZE['confine-bound'], 'STORE')
        self.prev_store_solv = os.path.join(self.dest, ORGANIZE['solvate-bound'], 'STORE')

        # [0.001, 0.01, 0.1, 0.2, 0.5, 1.0]
        # [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
        try:
            self.config = self.config['vba-bound']
        except KeyError:
            self.config = self.config['main']

        self.ll_list = list(map(lambda x: float(x), self.config.get('windows', '0.001, 0.01, 0.1, 0.2, 0.5, 1.0').split(', ')))
        self.int_meth = self.config.get('int_meth', 'TI')
        self.wham_path = self.config.get('WHAM_PATH', '')
        self.flucts = []
        self.dG = []

    def wham_write_timeseries(self, name, fluct):
        ll_list = [0] + self.ll_list
        for i in range(len(ll_list)):
            file_name = 'timeseries-' + name + '-' + str(ll_list[i])
            perc = 1
            list_length = int(len(fluct[i])*perc)
            f = open(os.path.join('WHAM', file_name), 'w')
            f.writelines(list(map(lambda x: str(x[0]) + ' ' + str(x[1]) + '\n', fluct[i][:list_length])))
            f.close()

    def wham_write_metadatafile(self, name, nb):
        ll_list = [0] + self.ll_list
        lines = []
        for i in range(len(ll_list)):
            lines.append('WHAM/timeseries-' + name + '-' + str(ll_list[i]) + ' ' + str(self.x0[nb]) + ' ' + str(float(self.kappa_max[nb])*ll_list[i]))
        file_name = 'metadatafile-' + name
        f = open(os.path.join('WHAM', file_name), 'w')
        f.writelines(list(map(lambda x: str(x) + '\n', lines)))
        f.close()

    def run(self):
        super(VbaBound, self).run()

        if not os.path.exists('STORE'):
            os.makedirs('STORE')

        # if not os.path.isfile('STORE/1.0.vba'):
        # Copy the topol file here
        shutil.copy(os.path.join(self.prev_store_solv, 'topol-complex-solv.top'), self.directory)

        filedata = self.create_plumed_tmpl()

        prev = os.path.join(self.prev_store, '6')
        for ll in self.ll_list:
            if not os.path.isfile('STORE/' + str(ll) + '.vba'):

                newdata = filedata.replace('KK1', str(self.kappa_max[0] * ll))
                newdata = newdata.replace('KK2', str(self.kappa_max[1] * ll))
                newdata = newdata.replace('KK3', str(self.kappa_max[2] * ll))
                newdata = newdata.replace('KK4', str(self.kappa_max[3] * ll))
                newdata = newdata.replace('KK5', str(self.kappa_max[4] * ll))
                newdata = newdata.replace('KK6', str(self.kappa_max[5] * ll))

                f = open('plumed.dat', 'w')
                f.write(newdata)
                f.close()

                subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, 'PRODUCTION.mdp') + ' -c ' + prev + '.gro -t ' + prev + '.cpt -p topol-complex-solv.top -o ' + str(ll) + '.tpr -maxwarn 2',
                                shell=True)
                subprocess.call('gmx mdrun -deffnm ' + str(ll) + ' -plumed plumed.dat -v',
                                shell=True)
                check_step('VBA_rest')
                shutil.move('VBA_rest', 'STORE/' + str(ll) + '.vba')
                self.store_files([str(ll) + '.gro', str(ll) + '.cpt'])
                prev = 'STORE/' + str(ll)
                os.remove('plumed.dat')
            else:
                prev = 'STORE/' + str(ll)
        clean_md_files()

        if self.int_meth == 'TI':
            # if not os.path.isfile('STORE/POS-ORIE.dhdl'):
            for ll in [0] + self.ll_list:
                fluct = str(ll)
                if ll == 0:
                    cols = [2, 3, 4, 5, 6, 7]
                    file = '../monitor-CVs/STORE/COLVAR-rest'
                else:
                    cols = [2, 4, 6, 8, 10, 12]
                    file = 'STORE/' + str(ll) + '.vba'
                for col in range(6):
                    fluct += ' ' + str(compute_fluct(self.x0[col], self.kappa_max[col], cols[col], file, dihe=self.dihe[col]))
                self.flucts.append(fluct)
            f = open('STORE/POS-ORIE.dhdl', 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', self.flucts)))
            f.close()

            if not self.flucts:
                with open('STORE/POS-ORIE.dhdl', 'r') as file:
                    for line in file:
                        self.flucts.append(line.rstrip('\n'))

            for i in range(2, 8):
                self.dG.append(compute_trapez(self.flucts, i))
            f = open('STORE/VBA_BND.dG', 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', self.dG)))
            f.close()

        if self.int_meth == 'WHAM':
            if not os.path.isdir("WHAM"):
                os.makedirs("WHAM")
            r_fluct = []
            tt_fluct = []
            phi_fluct = []
            TT_fluct = []
            PHI_fluct = []
            PSI_fluct = []
            for ll in [0] + self.ll_list:
                r_windows_fluct = []
                tt_windows_fluct = []
                phi_windows_fluct = []
                TT_windows_fluct = []
                PHI_windows_fluct = []
                PSI_windows_fluct = []
                if ll == 0:
                    file = '../monitor-CVs/STORE/COLVAR-rest'
                    with open(file, 'r') as f:
                        for line in f:
                            if not line.startswith('#'):
                                c = line.lstrip(' ').rstrip('\n').split(' ')
                                r_windows_fluct.append([float(c[0]), float(c[1])])
                                tt_windows_fluct.append([float(c[0]), float(c[2])])
                                phi_windows_fluct.append([float(c[0]), float(c[3])])
                                TT_windows_fluct.append([float(c[0]), float(c[4])])
                                PHI_windows_fluct.append([float(c[0]), float(c[5])])
                                PSI_windows_fluct.append([float(c[0]), float(c[6])])
                else:
                    file = 'STORE/' + str(ll) + '.vba'
                    with open(file, 'r') as f:
                        for line in f:
                            if not line.startswith('#'):
                                c = line.lstrip(' ').rstrip('\n').split(' ')
                                r_windows_fluct.append([float(c[0]), float(c[1])])
                                tt_windows_fluct.append([float(c[0]), float(c[3])])
                                phi_windows_fluct.append([float(c[0]), float(c[5])])
                                TT_windows_fluct.append([float(c[0]), float(c[7])])
                                PHI_windows_fluct.append([float(c[0]), float(c[9])])
                                PSI_windows_fluct.append([float(c[0]), float(c[11])])
                r_fluct.append(r_windows_fluct)
                tt_fluct.append(tt_windows_fluct)
                phi_fluct.append(phi_windows_fluct)
                TT_fluct.append(TT_windows_fluct)
                PHI_fluct.append(PHI_windows_fluct)
                PSI_fluct.append(PSI_windows_fluct)

            # r, tt, phi, TT, PHI, PSY
            for (name, fluct, nb) in [('r', r_fluct, 0), ('tt', tt_fluct, 1), ('phi', phi_fluct, 2), ('TT', TT_fluct, 3), ('PHI', PHI_fluct, 4), ('PSI', PSI_fluct, 5)]:
                self.wham_write_timeseries(name, fluct)
                self.wham_write_metadatafile(name, nb)
                command = os.path.join(self.wham_path, 'wham')
                mini = min([min([c[1] for c in fluct_windows]) for fluct_windows in fluct])
                maxi = max([max([c[1] for c in fluct_windows]) for fluct_windows in fluct])
                if name == 'r':
                    subprocess.call([command, str(mini), str(maxi), '60', '0.0000001', str(self.temp), '0', 'WHAM/metadatafile-' + name, 'wham-outfile-' + name, '600', '123456'])

                else:
                    # subprocess.call(command)
                    subprocess.call([command, 'Ppi', str(mini), str(maxi), '60', '0.0000001', str(self.temp), '0', 'WHAM/metadatafile-' + name, 'wham-outfile-' + name, '600', '123456'])

                output = open('wham-outfile-' + name, 'r')
                file_lines = output.readlines()
                output.close()
                dG = float(file_lines[-1].split()[1]) / 4.184
                self.dG.append(dG)

            f = open('STORE/VBA_BND.dG', 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', self.dG)))
            f.close()

        if not self.dG:
            with open('STORE/VBA_BND.dG', 'r') as file:
                for line in file:
                    self.dG.append(float(line.rstrip('\n')))

        return self.dG


def compute_sym_corr(sigma_l, sigma_p, sigma_pl, temp):
    kT = c.Boltzmann * c.N_A * temp / 1000  # kJ/mol
    return -kT * math.log(sigma_pl / (sigma_l * sigma_p)) / 4.184  # kcal/mol


class VbaUnbound(DDMClass):
    def __init__(self, config, kappa_max):
        super(VbaUnbound, self).__init__(config)
        self.directory = os.path.join(self.dest, ORGANIZE['vba-unbound'])
        self.prev_store = os.path.join(self.dest, ORGANIZE['vba-bound'], 'STORE')

        try:
            self.config = self.config['vba-unbound']
        except KeyError:
            self.config = self.config['main']

        self.sym_numbers = list(map(lambda x: float(x), self.config.get('symmetry_numbers', '1000, 1000, 1000').split(', ')))

        self.kappa = kappa_max
        self.dG = []
        self.sym_corr = []

    def run(self):
        super(VbaUnbound, self).run()

        if not os.path.exists('STORE'):
            os.makedirs('STORE')

        if not os.path.isfile('STORE/VBA_UB.dG'):
            for col in [2, 4, 8]:
                self.kappa.append(compute_mean(col, os.path.join(self.prev_store, '1.0.vba')))

            self.dG.append(compute_work(self.kappa, self.temp))
            f = open('STORE/VBA_UB.dG', 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', self.dG)))
            f.close()

        if not self.dG:
            with open('STORE/VBA_UB.dG', 'r') as file:
                for line in file:
                    self.dG.append(float(line.rstrip('\n')))

        if not os.path.isfile('STORE/sym_corr.dat'):
            if self.sym_numbers == [1000, 1000, 1000]:
                print("Vba unbound info : enable to compute symmetry correction, no symmetry number provided.")
                self.sym_corr = 0
            else:
                self.sym_corr.append(compute_sym_corr(1, 14, 1, self.temp))
            f = open('STORE/sym_corr.dat', 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', self.sym_corr)))
            f.close()

        if not self.sym_corr:
            with open('STORE/sym_corr.dat', 'r') as file:
                for line in file:
                    self.sym_corr.append(float(line.rstrip('\n')))

        return self.dG, self.sym_corr

