# -*- coding: utf-8 -*-
import os
import shutil
import subprocess

from .base import DDMClass, check_step, clean_md_files, compute_std, \
    compute_kf, compute_kf_plus, compute_fluct, compute_trapez, ORGANIZE


class Confine(DDMClass):
    def __init__(self, config, guest):
        super(Confine, self).__init__(config)

        self.guest = guest
        self.ll_list = [0.001, 0.01, 0.1, 0.2, 0.5, 1.0]

    def create_pdb_ref(self, input_file, output_file, renumbering=False):
        f = open(input_file, 'r')
        filedata = f.readlines()
        f.close()

        lines = []
        c = 1
        for line in filedata:
            if ' ' + self.guest.name + ' ' in line:
                str_c = str(c)
                while len(str_c) < 5:
                    str_c = ' ' + str_c

                line = line.rstrip('\n').replace('0.00 ', '1.00 ')
                if renumbering:
                    beg = line[:6]
                    end = line[12:]
                    line = beg + str_c + ' ' + end
                c += 1
                lines.append(line)

        f = open(output_file, 'w')
        f.writelines(list(map(lambda x: str(x) + '\n', lines)))
        f.close()

    def create_plumed_tmpl(self, type):
        f = open(os.path.join(self.static_dir, 'rmsd_' + type + '.tmpl'), 'r')
        filedata = f.read()
        f.close()

        newdata = filedata.replace('XXXXX', self.guest.beg)
        newdata = newdata.replace('YYYYY', str(int(self.guest.beg) + int(self.guest.nb_at) - 1))
        newdata = newdata.replace('RRRRR', self.guest.name)

        return newdata


class ConfineBound(Confine):
    def __init__(self, config, guest):
        super(ConfineBound, self).__init__(config, guest)
        self.directory = os.path.join(self.dest, ORGANIZE['confine-bound'])
        self.static_dir = os.path.join(self.static_dir, 'confine-bound')
        self.prev_store = os.path.join(self.dest, ORGANIZE['monitor-CVs'], 'STORE')
        self.prev_store_ref = os.path.join(self.dest, ORGANIZE['pick-reference'], 'STORE')
        self.prev_store_solv = os.path.join(self.dest, ORGANIZE['solvate-bound'], 'STORE')

        try:
            self.config = self.config['confine-bound']
        except KeyError:
            pass

        self.krms = ''
        self.krms_max = ''
        self.flucts = []
        self.dG = []

    def run(self):
        super(ConfineBound, self).run()

        # Extract coords of the ligand
        if not os.path.isfile('STORE/' + self.guest.name + '_ref.pdb'):
            self.create_pdb_ref(os.path.join(self.prev_store_ref, 'REFERENCE.pdb'), self.guest.name + '_ref.pdb')
            check_step(self.guest.name + '_ref.pdb')
            self.store_files([self.guest.name + '_ref.pdb'])

        # Monitor CV
        if not os.path.isfile('STORE/0.rms'):

            filedata = self.create_plumed_tmpl('anal')
            f = open(os.path.join(self.directory, 'plumed_rmsd_anal.inp'), 'w')
            f.write(filedata)
            f.close()

            subprocess.call('plumed driver --plumed plumed_rmsd_anal.inp --mf_xtc ' + os.path.join(self.prev_store_solv, 'prod.xtc') + ' --timestep 0.002 --trajectory-stride 2500',
                            shell=True)
            check_step('PLUMED.out')
            shutil.copy('PLUMED.out', 'STORE/0.rms')

        # Evaluate k_unbiased by QHA and kf
        if not os.path.isfile('STORE/file_max.krms'):
            # if not in config, compute.
            if not self.config.get('rms', False):
                std_c1 = compute_std(2, 'STORE/0.rms')

                # Compute Kf
                self.krms = compute_kf(std_c1, self.temp)
                self.krms_max = float(compute_kf_plus(self.krms))
            # if in config, store
            else:
                self.krms_max = float(self.config['rms'])

            f = open('STORE/file_max.krms', 'w')
            f.write(str(self.krms_max) + '\n')
            f.close()

        if self.krms_max == '':
            f = open('STORE/file_max.krms', 'r')
            kf = float(f.read())
            self.krms_max = kf
            f.close()

        if not os.path.isfile('STORE/6.gro'):
            # Copy the topol file here
            shutil.copy(os.path.join(self.prev_store_solv, 'topol-complex-solv.top'), self.directory)

            filedata = self.create_plumed_tmpl('bias')

            nn = 1
            prev = os.path.join(self.prev_store_solv, 'prod')
            for ll in self.ll_list :
                if not os.path.isfile('STORE/' + str(ll) + '.rms'):
                    kk = self.krms_max * ll

                    newdata = filedata.replace('KKKKK', str(kk))

                    f = open('plumed.dat', 'w')
                    f.write(newdata)
                    f.close()

                    subprocess.call('gmx grompp -f '+ os.path.join(self.static_dir, 'PRODUCTION.mdp') + ' -c ' + prev + '.gro -t ' + prev + '.cpt -p topol-complex-solv.top -o ' + str(nn) + '.tpr -maxwarn 2',
                                     shell=True)
                    subprocess.call('gmx_d mdrun -deffnm ' + str(nn) + ' -plumed plumed.dat -v',
                                    shell=True)
                    check_step('PLUMED-rmsd')
                    shutil.move('PLUMED-rmsd', 'STORE/' + str(ll) + '.rms')
                    shutil.move(str(nn) + '.gro', 'STORE/')
                    shutil.move(str(nn) + '.cpt', 'STORE/')
                    prev = 'STORE/' + str(nn)
                else:
                    prev = 'STORE/' + str(nn)
                nn += 1

            os.remove('plumed.dat')
            clean_md_files()

        if not os.path.isfile('STORE/RMS'):
            for ll in [0] + self.ll_list:
                self.flucts.append(str(ll) + ' ' + str(compute_fluct(0.0, self.krms_max, 2, "STORE/" + str(ll) + '.rms')))
            f = open('STORE/RMS', 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', self.flucts)))
            f.close()

        if not self.flucts:
            with open('STORE/RMS', 'r') as file:
                for line in file:
                    self.flucts.append(line.rstrip('\n'))

        if not os.path.isfile('STORE/CONF_BND.dG'):
            self.dG.append(compute_trapez(self.flucts, 2))
            f = open('STORE/CONF_BND.dG', 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', self.dG)))
            f.close()

        if not self.dG:
            with open('STORE/CONF_BND.dG', 'r') as file:
                for line in file:
                    self.dG.append(float(line.rstrip('\n')))

        return self.dG


class ConfineUnbound(Confine):
    def __init__(self, config, guest, krms_max):
        super(ConfineUnbound, self).__init__(config, guest)
        self.directory = os.path.join(self.dest, ORGANIZE['confine-unbound'])
        self.static_dir = os.path.join(self.static_dir, 'confine-unbound')
        self.prev_store = os.path.join(self.dest, ORGANIZE['monitor-CVs'], 'STORE')
        self.prev_store_ref = os.path.join(self.dest, ORGANIZE['pick-reference'], 'STORE')
        self.prev_store_solv = os.path.join(self.dest, ORGANIZE['solvate-unbound'], 'STORE')

        self.krms_max = krms_max
        self.flucts = []
        self.dG = []

    def run(self):
        super(ConfineUnbound, self).run()

        if not os.path.isfile('STORE/' + self.guest.name + '_ref.pdb'):
            self.create_pdb_ref(os.path.join(self.prev_store_ref, 'REFERENCE.pdb'), self.guest.name + '_ref.pdb', renumbering=True)
            check_step(self.guest.name + '_ref.pdb')
            self.store_files([self.guest.name + '_ref.pdb'])

        # Monitor CV
        if not os.path.isfile('STORE/0.rms'):
            self.guest.beg = '1'
            filedata = self.create_plumed_tmpl('anal')
            f = open(os.path.join(self.directory, 'plumed_rmsd_anal.inp'), 'w')
            f.write(filedata)
            f.close()

            subprocess.call('plumed driver --plumed plumed_rmsd_anal.inp --mf_xtc ' + os.path.join(self.prev_store_solv, 'prod.xtc') + ' --timestep 0.002 --trajectory-stride 2500',
                            shell=True)
            check_step('PLUMED.out')
            shutil.move('PLUMED.out', 'STORE/0.rms')

        if not os.path.isfile('STORE/1.0.rms'):
            # Copy the topol file here
            shutil.copy(os.path.join(self.prev_store_solv, 'topol-ligand-solv.top'), self.directory)

            self.guest.beg = '1'
            filedata = self.create_plumed_tmpl('bias')

            nn = 1
            prev = os.path.join(self.prev_store_solv, 'prod')
            for ll in self.ll_list:
                if not os.path.isfile('STORE/' + str(ll) + '.rms'):
                    kk = self.krms_max * ll

                    newdata = filedata.replace('KKKKK', str(kk))

                    f = open('plumed.dat', 'w')
                    f.write(newdata)
                    f.close()

                    subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, 'PRODUCTION.mdp') + ' -c ' + prev + '.gro -t ' + prev + '.cpt -p topol-ligand-solv.top -o ' + str(nn) + '.tpr -maxwarn 2',
                                     shell=True)
                    subprocess.call('gmx_d mdrun -deffnm ' + str(nn) + ' -plumed plumed.dat -v',
                                    shell=True)
                    check_step('PLUMED-rmsd')
                    shutil.move('PLUMED-rmsd', 'STORE/' + str(ll) + '.rms')
                    shutil.move(str(nn) + '.gro', 'STORE/')
                    shutil.move(str(nn) + '.cpt', 'STORE/')
                    prev = 'STORE/' + str(nn)
                else:
                    prev = 'STORE/' + str(nn)
                nn += 1

            os.remove('plumed.dat')
            clean_md_files()

        if not os.path.isfile('STORE/RMS'):
            for ll in [0] + self.ll_list:
                self.flucts.append(str(ll) + ' ' + str(compute_fluct(0.0, self.krms_max, 2, "STORE/" + str(ll) + '.rms')))
            f = open('STORE/RMS', 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', self.flucts)))
            f.close()

        if not self.flucts:
            with open('STORE/RMS', 'r') as file:
                for line in file:
                    self.flucts.append(line.rstrip('\n'))

        if not os.path.isfile('STORE/CONF_BND.dG'):
            self.dG.append(compute_trapez(self.flucts, 2))
            f = open('STORE/CONF_BND.dG', 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', self.dG)))
            f.close()

        if not self.dG:
            with open('STORE/CONF_BND.dG', 'r') as file:
                for line in file:
                    self.dG.append(float(line.rstrip('\n')))

        return self.dG
