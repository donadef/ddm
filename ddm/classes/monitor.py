# -*- coding: utf-8 -*-
import os
import subprocess

from .base import DDMClass, check_step, compute_std, compute_kf, compute_kf_plus


class MonitorCVs(DDMClass):
    def __init__(self, config, complex):
        super(MonitorCVs, self).__init__(config, complex)

        self.prev_store = os.path.join(self.dest, '03-pick-reference/STORE')
        self.prev_store_solv = os.path.join(self.dest, '01-solvate-bound/STORE')
        self.directory = os.path.join(self.dest, '04-monitor-CVs')

        self.x0 = []
        self.kappa = []
        self.krms = []

    def run(self):
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

        os.chdir(self.directory)

        # Monitor POS/ORIE of ligand in REFERENCE
        if not os.path.isfile('STORE/file.x0'):
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
        if not os.path.isfile('STORE/file.kappa'):
            if not os.path.exists('STORE'):
                os.makedirs('STORE')

            subprocess.call('plumed driver --plumed ' + os.path.join(self.prev_store, 'vba.dat') + ' --mf_xtc ' + os.path.join(self.prev_store_solv, 'prod.xtc') + ' --timestep 0.002',
                            shell=True)

            check_step('COLVAR-rest')

            # Plot distributions
            # nb_bin = subprocess.check_output("wc COLVAR-rest | awk '{print sqrt($1)}'", shell=True)
            # subprocess.call('awk -f ' + os.path.join(self.awk_dir, 'histo.awk') + ' -v col=2 -v NBIN=' + nb_bin + ' COLVAR-rest > rr.hh', shell=True)
            # subprocess.call('awk -f ' + os.path.join(self.awk_dir, 'histo.awk') + ' -v col=3 -v NBIN=' + nb_bin + ' COLVAR-rest > tt.hh', shell=True)
            # subprocess.call('awk -f ' + os.path.join(self.awk_dir, 'histo.awk') + ' -v col=4 -v NBIN=' + nb_bin + ' COLVAR-rest > phi.hh', shell=True)
            # subprocess.call('awk -f ' + os.path.join(self.awk_dir, 'histo.awk') + ' -v col=5 -v NBIN=' + nb_bin + ' COLVAR-rest > TT.hh', shell=True)
            # subprocess.call('awk -f ' + os.path.join(self.awk_dir, 'histo.awk') + ' -v col=6 -v NBIN=' + nb_bin + ' COLVAR-rest > PHI.hh', shell=True)
            # subprocess.call('awk -f ' + os.path.join(self.awk_dir, 'histo.awk') + ' -v col=7 -v NBIN=' + nb_bin + ' COLVAR-rest > PSI.hh', shell=True)

            std_cvs = []
            for col in range(2, 8):
                std_cvs.append(compute_std(col, 'COLVAR-rest'))

            # Compute Kfs
            self.kappa = list(map(compute_kf, std_cvs))
            f = open('STORE/file.kappa', 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', self.kappa)))
            f.close()

            os.remove('COLVAR-rest')
            check_step('STORE/file.kappa')

        if not self.kappa:
            with open('STORE/file.kappa', 'r') as file:
                for line in file:
                    self.kappa.append(float(line.rstrip('\n')))

        if not os.path.isfile('STORE/file.krms'):
            self.krms = list(map(compute_kf_plus, self.kappa))
            f = open('STORE/file.krms', 'w')
            f.writelines(list(map(lambda x: str(x) + '\n', self.krms)))
            f.close()

        if not self.krms:
            with open('STORE/file.krms', 'r') as file:
                for line in file:
                    self.krms.append(float(line.rstrip('\n')))

        self.store_files()
