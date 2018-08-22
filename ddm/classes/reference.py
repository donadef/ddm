# -*- coding: utf-8 -*-
import os
import subprocess

from .base import DDMClass


class PickReference(DDMClass):
    def __init__(self, config, complex):
        super(PickReference, self).__init__(config, complex)

        self.prev_store = os.path.join(self.dest, '00-modeling/STORE')
        self.prev_store_solv = os.path.join(self.dest, '01-solvate-bound/STORE')
        self.directory = os.path.join(self.dest, '03-pick-reference')
        self.static_dir = os.path.join(self.static_dir, '03-pick-reference')

    def run(self):
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

        os.chdir(self.directory)

        # Extract last 10 ns (recentered)
        if not os.path.isfile('last10ns.xtc'):
            subprocess.call('echo "Other Other" | gmx trjconv  -f ' + os.path.join(self.prev_store_solv, 'prod.xtc') + ' -b 10000 -o last10ns.xtc -pbc cluster -s ' + os.path.join(self.prev_store_solv, 'prod.tpr'),
                            shell=True)
        self.check_step('last10ns.xtc')

        # Cluster trj based on RMSD (cutoff 0.03 nm)
        if not os.path.isfile('clusters.pdb'):
            subprocess.call('echo "Other Other" | gmx cluster -f last10ns.xtc -s prod.tpr -dist -cutoff 0.03 -cl -method jarvis-patrick',
                            shell=True)
        self.check_step('clusters.pdb')

        if not os.path.isfile('most_populated.pdb'):
            nb_row = subprocess.check_output("grep ATOM " + os.path.join(self.prev_store, 'complex_mini.pdb') + " | tail -1 | awk '{print $2+2}'",
                                             shell=True).decode("utf-8").rstrip('\n')
            subprocess.call('head -' + nb_row + ' clusters.pdb > most_populated.pdb')
        self.check_step('most_populated.pdb')

        nb_step = ''
        if not os.path.isfile('relaxed.trr'):
            subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, 'MINI.mdp') + ' -c most_populated.pdb -p ' + os.path.join(self.prev_store, 'topol-complex.top') + ' -o relaxed -maxwarn 2',
                            shell=True)
            subprocess.call('gmx mdrun -v -deffnm relaxed > TMP 2>&1',
                            shell=True)
        self.check_step('relaxed.trr')

        nb_step = subprocess.check_output("grep 'Step' TMP | tail -1 | awk '{print $2}' | sed s/','//g",
                                          shell=True).decode("utf-8").rstrip('\n')
        subprocess.call('echo "System" | gmx trjconv -f relaxed.trr -s relaxed.tpr -o REFERENCE.pdb -b ' + nb_step,
                        shell=True)

        self.clean_tmp()

        if not os.path.isfile('vba.dat'):
            f = open(os.path.join(self.static_dir, 'vba.tmpl'), 'r')
            filedata = f.read()
            f.close()

            newdata = filedata.replace('P1', self.config['03-pick-reference']['P1'])
            newdata = newdata.replace('P2', self.config['03-pick-reference']['P2'])
            newdata = newdata.replace('P3', self.config['03-pick-reference']['P3'])
            newdata = newdata.replace('L1', self.config['03-pick-reference']['L1'])
            newdata = newdata.replace('L2', self.config['03-pick-reference']['L2'])
            newdata = newdata.replace('L3', self.config['03-pick-reference']['L3'])

            f = open(os.path.join(self.directory, 'STORE/vba.dat'), 'w')
            f.write(newdata)
            f.close()

        files_to_store = ['REFERENCE.pdb']
        self.store_files(files_to_store)
