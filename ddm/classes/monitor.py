# -*- coding: utf-8 -*-
import os
import subprocess

from .base import DDMClass


class MonitorCVs(DDMClass):
    def __init__(self, config, complex):
        super(MonitorCVs, self).__init__(config, complex)

        self.prev_store = os.path.join(self.dest, '03-pick-reference/STORE')
        self.prev_store_solv = os.path.join(self.dest, '01-solvate-bound/STORE')
        self.directory = os.path.join(self.dest, '04-monitor-CVs')

    def run(self):
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

        os.chdir(self.directory)

        # Monitor POS/ORIE of ligand in REFERENCE
        if not os.path.isfile('file.x0'):
            subprocess.call('plumed driver --plumed ' + os.path.join(self.prev_store, 'vba.dat') + ' --mf_pdb REFERENCE.pdb',
                            shell=True)
            # TODO : Do this using python, not awk !! Read lines and write new file.
            subprocess.call("awk 'substr($0,1,1)!=\"#\"{for(ii=2;ii<=NF;ii++){print $ii}}' COLVAR-rest > file.x0 ",
                            shell=True)
        self.check_step('file.x0')
        os.remove('COLVAR-rest')

        # Monitor POS and ORIE of ligand in unbiased MD
        if not os.path.isfile('file.kappa'):
            subprocess.call('plumed driver --plumed ' + os.path.join(self.prev_store, 'vba.dat') + ' --mf_xtc prod.xtc --timestep 0.002',
                            shell=True)

            self.check_step('COLVAR-rest')

            # Plot distributions
            # nb_bin = subprocess.check_output("wc COLVAR-rest | awk '{print sqrt($1)}'", shell=True)
            # subprocess.call('awk -f ' + os.path.join(self.awk_dir, 'histo.awk') + ' -v col=2 -v NBIN=' + nb_bin + ' COLVAR-rest > rr.hh', shell=True)
            # subprocess.call('awk -f ' + os.path.join(self.awk_dir, 'histo.awk') + ' -v col=3 -v NBIN=' + nb_bin + ' COLVAR-rest > tt.hh', shell=True)
            # subprocess.call('awk -f ' + os.path.join(self.awk_dir, 'histo.awk') + ' -v col=4 -v NBIN=' + nb_bin + ' COLVAR-rest > phi.hh', shell=True)
            # subprocess.call('awk -f ' + os.path.join(self.awk_dir, 'histo.awk') + ' -v col=5 -v NBIN=' + nb_bin + ' COLVAR-rest > TT.hh', shell=True)
            # subprocess.call('awk -f ' + os.path.join(self.awk_dir, 'histo.awk') + ' -v col=6 -v NBIN=' + nb_bin + ' COLVAR-rest > PHI.hh', shell=True)
            # subprocess.call('awk -f ' + os.path.join(self.awk_dir, 'histo.awk') + ' -v col=7 -v NBIN=' + nb_bin + ' COLVAR-rest > PSI.hh', shell=True)

            lcv1 = []
            lcv2 = []
            lcv3 = []
            lcv4 = []
            lcv5 = []
            lcv6 = []
            with open('COLVAR-rest', 'r') as colvar_file:
                for line in colvar_file:
                    if not line.startswith('#'):
                        trash, cv1, cv2, cv3, cv4, cv5, cv6 = line.lstrip(' ').rstrip('\n').split(' ')
                        print(cv1, cv2, cv3, cv4, cv5, cv6)
