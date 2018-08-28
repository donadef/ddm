# -*- coding: utf-8 -*-
import os
import shutil
import subprocess

from .base import DDMClass, check_step, clean_md_files, compute_std, compute_fluct, compute_trapez, compute_work


class Vba(DDMClass):
    def __init__(self, config, guest, x0, kappa_max, krms_max):
        super(Vba, self).__init__(config)

        self.guest = guest
        self.x0 = x0
        self.kappa_max = kappa_max
        self.krms_max = krms_max

        self.dihe = [False, False, True, False, True, True]

    def create_plumed_tmpl(self):
        f = open(os.path.join(self.static_dir, 'vba_bias.tmpl'), 'r')
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


class VbaBound(Vba):
    def __init__(self, config, guest, x0, kappa_max, krms_max):
        super(VbaBound, self).__init__(config, guest, x0, kappa_max, krms_max)
        self.directory = os.path.join(self.dest, '06-vba-bound')
        self.static_dir = os.path.join(self.static_dir, '06-vba-bound')
        self.prev_store = os.path.join(self.dest, '05-confine-bound/STORE')
        self.prev_store_solv = os.path.join(self.dest, '01-solvate-bound/STORE')

        self.flucts = []
        self.dG = []

    def run(self):
        super(VbaBound, self).run()

        if not os.path.isfile('STORE/1.0.vba'):
            if not os.path.exists('STORE'):
                os.makedirs('STORE')

            # Copy the topol file here
            shutil.copy(os.path.join(self.prev_store_solv, 'topol-complex-solv.top'), self.directory)

            filedata = self.create_plumed_tmpl()

            nn = 1
            prev = os.path.join(self.prev_store, '6')
            for ll in [0.001, 0.01, 0.1, 0.2, 0.5, 1.0]:
                if not os.path.isfile('STORE' + str(ll) + '.vba'):

                    newdata = filedata.replace('KK1', str(self.kappa_max[0] * ll))
                    newdata = newdata.replace('KK2', str(self.kappa_max[1] * ll))
                    newdata = newdata.replace('KK3', str(self.kappa_max[2] * ll))
                    newdata = newdata.replace('KK4', str(self.kappa_max[3] * ll))
                    newdata = newdata.replace('KK5', str(self.kappa_max[4] * ll))
                    newdata = newdata.replace('KK6', str(self.kappa_max[5] * ll))

                    f = open('plumed.dat', 'w')
                    f.write(newdata)
                    f.close()

                    subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, 'PRODUCTION.mdp') + ' -c ' + prev + '.gro -t ' + prev + '.cpt -p topol-complex-solv.top -o ' + str(nn) + '.tpr -maxwarn 2',
                                    shell=True)
                    subprocess.call('gmx mdrun -deffnm ' + str(nn) + ' -plumed plumed.dat -v',
                                    shell=True)
                    check_step('VBA_rest')
                    shutil.move('VBA_rest', 'STORE/' + str(ll) + '.vba')
                    prev = str(nn)
                else:
                    prev = 'STORE/' + str(nn)
                nn += 1

            os.remove('plumed.dat')
        clean_md_files()

        if not os.path.isfile('STORE/POS-ORIE.dhdl'):
            for ll in [0, 0.001, 0.01, 0.1, 0.2, 0.5, 1.0]:
                fluct = str(ll)
                if ll == 0:
                    cols = [2, 3, 4, 5, 6, 7]
                    file = '../04-monitor-CVs/STORE/COLVAR-rest'
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

        if not os.path.isfile('STORE/VBA_BND.dG'):
            for i in range(2, 8):
                self.dG.append(compute_trapez(self.flucts, i))
            f = open('STORE/VBA_BND.dG', 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', self.dG)))
            f.close()

        if not self.dG:
            with open('STORE/VBA_BND.dG', 'r') as file:
                for line in file:
                    self.dG.append(line.rstrip('\n'))

        return self.dG


class VbaUnbound(DDMClass):
    def __init__(self, config, kappa_max):
        super(VbaUnbound, self).__init__(config)
        self.directory = os.path.join(self.dest, '07-work-unbound')
        self.prev_store = os.path.join(self.dest, '06-vba-bound/STORE')

        self.kappa_max = kappa_max
        self.dG = []

    def run(self):
        super(VbaUnbound, self).run()

        for col in [2, 3, 4]:
            self.kappa_max.append(compute_std(col, os.path.join(self.prev_store, '1.0.vba')))
        print(self.kappa_max)

        if not os.path.isfile('STORE/VBA_UB.dG'):
            if not os.path.exists('STORE'):
                os.makedirs('STORE')
            self.dG = compute_work(self.kappa_max)
            f = open('STORE/VBA_UB.dG', 'w')
            f.write(str(self.dG) + '\n')
            f.close()

        if self.dG == '':
            with open('STORE/VBA_UB.dG', 'r') as file:
                for line in file:
                    self.dG.append(line.rstrip('\n'))

        return self.dG
