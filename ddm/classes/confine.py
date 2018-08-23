# -*- coding: utf-8 -*-
import os
import shutil
import subprocess
import numpy as np

from .base import DDMClass, check_step


class Confine(DDMClass):
    def __init__(self, config, complex):
        super(Confine, self).__init__(config, complex)

    def run(self):
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

        os.chdir(self.directory)

    def create_plumed_tmpl(self, type, beg, end):
        f = open(os.path.join(self.static_dir, 'rmsd_' + type + '.tmpl'), 'r')
        filedata = f.read()
        f.close()

        newdata = filedata.replace('XXXXX', str(beg))
        newdata = newdata.replace('YYYYY', str(end))
        newdata = newdata.replace('RRRRR', self.guest)

        f = open(os.path.join(self.directory, 'plumed_rmsd_' + type + '.inp'), 'w')
        f.write(newdata)
        f.close()


class ConfineBound(Confine):
    def __init__(self, config, complex):
        super(ConfineBound, self).__init__(config, complex)
        self.directory = os.path.join(self.dest, '05-confine-bound')
        self.static_dir = os.path.join(self.static_dir, '05-confine-bound')
        self.prev_store = os.path.join(self.dest, '04-monitor-CVs/STORE')
        self.prev_store_ref = os.path.join(self.dest, '03-pick-reference/STORE')
        self.prev_store_solv = os.path.join(self.dest, '01-solvate-bound/STORE')

    def run(self):
        super(ConfineBound, self).run()

        # Extract coords of the ligand
        if not os.path.isfile(self.guest + '_ref.pdb'):
            subprocess.call("sed s/'0.00 '/'1.00 '/g " + os.path.join(self.prev_store_ref, 'REFERENCE.pdb') + " | grep " + self.guest + " > " + self.guest + "_ref.pdb",
                            shell=True)
            check_step(self.guest + '_ref.pdb')

        # Create the plumed template for rmsd
        if not os.path.isfile('plumed_rmsd_anal.inp') or not os.path.isfile('plumed_rmsd_bias.inp'):
            beg = subprocess.check_output("grep 'ATOM' " + self.guest + "_ref.pdb | head -1 | awk '{print $2}'",
                                          shell=True).decode("utf-8").rstrip('\n')
            end = subprocess.check_output("tail -1 " + self.guest + "_ref.pdb | awk '{print $2}'",
                                          shell=True).decode("utf-8").rstrip('\n')
            if not os.path.isfile('plumed_rmsd_anal.inp'):
                self.create_plumed_tmpl('anal', beg, end)
                check_step('plumed_rmsd_anal.inp')
            if not os.path.isfile('plumed_rmsd_bias.inp'):
                self.create_plumed_tmpl('bias', beg, end)
                check_step('plumed_rmsd_bias.inp')

        # Monitor CV
        path_to_plumed_out = ''
        if not os.path.isfile('PLUMED.out'):
            if not os.path.isfile('STORE/PLUMED.out'):
                subprocess.call('plumed driver --plumed plumed_rmsd_anal.inp --mf_xtc ' + os.path.join(self.prev_store_solv, 'prod.xtc') + ' --timestep 0.002 --trajectory-stride 2500',
                                shell=True)
                check_step('PLUMED.out')
                self.files_to_store.append('PLUMED.out')
            else:
                path_to_plumed_out = 'STORE/'

        self.store_files()

        kf = '1'
        # Evaluate k_unbiased by QHA and kf
        if not os.path.isfile('STORE/file.kappa') or not os.path.isfile('STORE/file.krms'):
            lc1 = []
            with open(path_to_plumed_out + 'PLUMED.out', 'r') as plumed_file:
                for line in plumed_file:
                    if not line.startswith('#'):
                        trash, c1 = line.lstrip(' ').rstrip('\n').split(' ')
                        lc1.append(float(c1))
            a = np.array(lc1)
            # Compute the standard deviation
            std_c1 = np.std(a)

            # Compute Kf
            kf_u = (8.314 * 298 / 1000) / std_c1 ** 2
            f = open('STORE/file.kappa', 'w')
            f.write(str(kf_u) + '\n')
            f.close()
            check_step('STORE/file.kappa')

            kf += '0' * len(str(kf_u).split('.')[0])
            f = open('STORE/file.krms', 'w')
            f.write(str(kf) + '\n')
            f.close()

        if kf == '1':
            f = open('STORE/file.krms', 'r')
            kf = f.read()
            f.close()

        if not os.path.isfile('STORE/6.gro'):
            # Copy the topol file here
            shutil.copy(os.path.join(self.prev_store_solv, 'topol-complex-solv.top'), self.directory)
            nn = 1
            prev = os.path.join(self.prev_store_solv, 'prod')
            for ll in [0.001, 0.01, 0.1, 0.2, 0.5, 1.0]:
                if not os.path.isfile('STORE' + str(ll) + '.rms'):
                    kk = float(kf) * ll

                    # Modify the plumed_rmsd_bias.inp file
                    f = open('plumed_rmsd_bias.inp', 'r')
                    filedata = f.read()
                    f.close()

                    newdata = filedata.replace('KKKKK', str(kk))

                    f = open('file.dat', 'w')
                    f.write(newdata)
                    f.close()

                    subprocess.call('gmx grompp -f '+ os.path.join(self.static_dir, 'PRODUCTION.mdp') + ' -c ' + prev + '.gro -t ' + prev + '.cpt -p topol-complex-solv.top -o ' + str(nn) + '.tpr -maxwarn 2',
                                     shell=True)
                    subprocess.call('gmx_d mdrun -deffnm ' + str(nn) + ' -plumed file.dat -v',
                                    shell=True)
                    check_step('PLUMED-rmsd')
                    shutil.move('PLUMED-rmsd', 'STORE/' + str(ll) + '.rms')
                prev = str(nn)
                nn += 1

            self.files_to_store = ['6.gro', '6.cpt']
            self.store_files()



