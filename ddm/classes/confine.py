# -*- coding: utf-8 -*-
import os
import shutil
import subprocess

from .base import DDMClass, check_step, clean_md_files, compute_std, compute_kf, compute_kf_plus, compute_fluct, compute_trapez, compute_work


class Confine(DDMClass):
    def __init__(self, config, guest):
        super(Confine, self).__init__(config)

        self.guest = guest

    def create_plumed_tmpl(self, type):
        f = open(os.path.join(self.static_dir, 'rmsd_' + type + '.tmpl'), 'r')
        filedata = f.read()
        f.close()

        newdata = filedata.replace('XXXXX', self.guest.beg)
        newdata = newdata.replace('YYYYY', self.guest.end)
        newdata = newdata.replace('RRRRR', self.guest.name)

        return newdata


class ConfineBound(Confine):
    def __init__(self, config, guest):
        super(ConfineBound, self).__init__(config, guest)
        self.directory = os.path.join(self.dest, '05-confine-bound')
        self.static_dir = os.path.join(self.static_dir, '05-confine-bound')
        self.prev_store = os.path.join(self.dest, '04-monitor-CVs/STORE')
        self.prev_store_ref = os.path.join(self.dest, '03-pick-reference/STORE')
        self.prev_store_solv = os.path.join(self.dest, '01-solvate-bound/STORE')

        self.krms = ''
        self.krms_max = ''
        self.flucts = []
        self.dG = []

    def run(self):
        super(ConfineBound, self).run()

        # Extract coords of the ligand
        if not os.path.isfile('STORE/' + self.guest.name + '_ref.pdb'):
            # TODO: use python !
            subprocess.call("sed s/'0.00 '/'1.00 '/g " + os.path.join(self.prev_store_ref, 'REFERENCE.pdb') + " | grep " + self.guest.name + " > " + self.guest.name + "_ref.pdb",
                            shell=True)
            check_step(self.guest.name + '_ref.pdb')
            self.files_to_store = [self.guest.name + '_ref.pdb']
            self.store_files()

        # Monitor CV
        path_to_plumed_out = ''
        if not os.path.isfile('PLUMED.out'):
            if not os.path.isfile('STORE/PLUMED.out'):

                filedata = self.create_plumed_tmpl('anal')
                f = open(os.path.join(self.directory, 'plumed_rmsd_anal.inp'), 'w')
                f.write(filedata)
                f.close()

                subprocess.call('plumed driver --plumed plumed_rmsd_anal.inp --mf_xtc ' + os.path.join(self.prev_store_solv, 'prod.xtc') + ' --timestep 0.002 --trajectory-stride 2500',
                                shell=True)
                check_step('PLUMED.out')
                shutil.copy('PLUMED.out', 'STORE/0.rms')
            else:
                path_to_plumed_out = 'STORE/'

        # Evaluate k_unbiased by QHA and kf
        if not os.path.isfile('STORE/file.kappa'):
            std_c1 = compute_std(2, path_to_plumed_out + 'PLUMED.out')

            # Compute Kf
            self.krms = compute_kf(std_c1)
            f = open('STORE/file.krms', 'w')
            f.write(str(self.krms) + '\n')
            f.close()
            check_step('STORE/file.krms')

        if self.krms == '':
            f = open('STORE/file.krms', 'r')
            self.krms = float(f.read())
            f.close()

        if not os.path.isfile('STORE/file_max.krms'):
            self.krms_max = float(compute_kf_plus(self.krms))
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
            for ll in [0.001, 0.01, 0.1, 0.2, 0.5, 1.0]:
                if not os.path.isfile('STORE' + str(ll) + '.rms'):
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
                    prev = str(nn)
                else:
                    prev = 'STORE/' + str(nn)
                nn += 1

            self.files_to_store = ['6.gro', '6.cpt']
            self.store_files()

            os.remove('plumed.dat')
            clean_md_files()

        if not os.path.isfile('STORE/RMS'):
            for ll in [0, 0.001, 0.01, 0.1, 0.2, 0.5, 1.0]:
                self.flucts.append(str(ll) + ' ' + str(compute_fluct(0.0, self.krms, 2, "STORE/" + str(ll) + '.rms')))
            f = open('STORE/RMS', 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', self.flucts)))
            f.close()

        if not self.flucts:
            with open('STORE/RMS', 'r') as file:
                for line in file:
                    self.flucts.append(line.rstrip('\n'))

        if not os.path.isfile('STORE/CONF_BND.dG'):
            self.dG = compute_trapez(self.flucts, 2)
            f = open('STORE/CONF_BND.dG', 'w')
            f.write(str(self.dG) + '\n')
            f.close()

        if not self.dG:
            with open('STORE/CONF_BND.dG', 'r') as file:
                for line in file:
                    self.dG.append(line.rstrip('\n'))

        return self.dG

