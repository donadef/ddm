# -*- coding: utf-8 -*-
import os
import subprocess

from .base import DDMClass, check_step, compute_std, compute_kf, compute_kf_plus, ORGANIZE


class MonitorCVs(DDMClass):
    def __init__(self, config):
        super(MonitorCVs, self).__init__(config)

        self.prev_store = os.path.join(self.dest, ORGANIZE['pick-reference'], 'STORE')
        self.prev_store_solv = os.path.join(self.dest, ORGANIZE['solvate-bound'], 'STORE')
        self.directory = os.path.join(self.dest, ORGANIZE['monitor-CVs'])
        self.config = self.config['monitor-CVs']

        self.x0 = []
        self.kappa = []
        self.kappa_max = []

    def run(self):
        super(MonitorCVs, self).run()

        # Monitor POS/ORIE of ligand in REFERENCE
        if not os.path.isfile('STORE/file.x0'):
            if not os.path.exists('STORE'):
                os.makedirs('STORE')

            subprocess.call('plumed driver --plumed ' + os.path.join(self.prev_store, 'vba.dat') + ' --mf_pdb ' + os.path.join(self.prev_store, 'REFERENCE.pdb'),
                            shell=True)

            check_step('COLVAR-rest')

            with open('COLVAR-rest', 'r') as file:
                for line in file:
                    if not line.startswith('#'):
                        trash, c1, c2, c3, c4, c5, c6 = line.lstrip(' ').rstrip('\n').split(' ')
                        self.x0 = [c1, c2, c3, c4, c5, c6]
            f = open('STORE/file.x0', 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', self.x0)))
            f.close()

            os.remove('COLVAR-rest')
            check_step('STORE/file.x0')

        if not self.x0:
            with open('STORE/file.x0', 'r') as file:
                for line in file:
                    self.x0.append(float(line.rstrip('\n')))


        # Monitor POS and ORIE of ligand in unbiased MD
        if not os.path.isfile('STORE/file_max.kappa'):
            # if the over-estimated restrains are not in the config file, compute them
            if not self.config['rr'] or not self.config['tt'] or not self.config['phi'] or not self.config['TT'] or not self.config['PHI'] or not self.config['PSI']:
                subprocess.call('plumed driver --plumed ' + os.path.join(self.prev_store, 'vba.dat') + ' --mf_xtc ' + os.path.join(self.prev_store_solv, 'prod.xtc') + ' --timestep 0.002',
                                shell=True)

                check_step('COLVAR-rest')
                self.files_to_store = ['COLVAR-rest']
                self.store_files()

                std_cvs = []
                for col in range(2, 8):
                    std_cvs.append(compute_std(col, 'COLVAR-rest'))

                # Compute Kfs
                self.kappa = list(map(compute_kf, std_cvs))

                os.remove('COLVAR-rest')
                self.kappa_max = list(map(compute_kf_plus, self.kappa))
            # if the over-estimated restrains are in the config file, just save them in self.kappa_max
            else:
                self.kappa_max = [self.config['rr'], self.config['tt'], self.config['phi'],
                                  self.config['TT_'], self.config['PHI_'], self.config['PSI']]
            if not os.path.exists('STORE'):
                os.makedirs('STORE')
            f = open('STORE/file_max.kappa', 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', self.kappa_max)))
            f.close()

        if not self.kappa_max:
            with open('STORE/file_max.kappa', 'r') as file:
                for line in file:
                    self.kappa_max.append(float(line.rstrip('\n')))
