# -*- coding: utf-8 -*-
import os
import shutil
import subprocess
from tools import cgenff_charmm2gmx


class Modeling:
    def __init__(self, host, guest, complex, ipdb, dest, static_dir):
        self.host = host
        self.guest = guest
        self.complex = complex
        self.ipdb = ipdb
        self.dest = dest
        self.static_dir = static_dir

        self.directory = os.path.join(self.dest, '00-modeling')

    def run(self):
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

        # Copy the complex pdb in the folder
        shutil.copy(self.complex, self.directory)

        os.chdir(self.directory)

        # Guest
        self.prepare_guest()
        self.minimize_guest()

        # Host
        self.prepare_host()

        # Complex
        self.minimize_complex()

        files_to_store = [self.guest + '.top', self.guest + '.prm', self.guest + '.itp', self.guest + '_ini.pdb', 'topol-ligand.top',
                          self.host + '.top', self.host + '.prm', self.host + '.itp', self.host + '_ini.pdb',
                          'complex_mini.pdb', 'topol-complex.top']
        store_dir = os.path.join(self.directory, 'STORE')
        if not os.path.exists(store_dir):
            os.makedirs(store_dir)

        for file in files_to_store:
            try:
                shutil.copy(file, store_dir)
            except FileNotFoundError:
                print('File ' + file + ' has not been created.')
                exit()

    def prepare_guest(self):
        # PDB file for the guest, extract from the complex given as input.
        if not os.path.isfile(os.path.join(self.directory, self.guest + '.pdb')):
            subprocess.call('grep " ' + self.guest + ' " ' + self.ipdb + '.pdb' + ' > ' + os.path.join(self.directory, self.guest + '.pdb'),
                            shell=True)
        if not os.path.isfile(os.path.join(self.directory, self.guest + '.top')):
            self.generate_param_files(self.guest)

    def minimize_guest(self):
        if not os.path.isfile(os.path.join(self.directory, 'topol-ligand.top')):
            # create the topol-ligand.top file
            f = open(os.path.join(self.static_dir, '00-modeling/ligand.top'), 'r')
            filedata = f.read()
            f.close()

            newdata = filedata.replace("XXXXX", self.guest)

            f = open(os.path.join(self.directory, 'topol-ligand.top'), 'w')
            f.write(newdata)
            f.close()

        # First minimization step, steepest descent method.
        if not os.path.isfile(os.path.join(self.directory, 'mini1.trr')):
            self.minimize1(self.guest + '_ini', 'topol-ligand')

        # Second minimization step, conjugate gradient method.
        if not os.path.isfile(os.path.join(self.directory, 'mini2.trr')):
            self.minimize2(self.guest + '_ini', 'topol-ligand')

        # Extract coordinates
        if not os.path.isfile(os.path.join(self.directory, 'ligand_mini.pdb')):
            if not os.path.isfile(os.path.join(self.directory, 'TMP')):
                self.minimize2(self.guest + '_ini', 'topol-ligand')
            self.extract_coordinates('ligand_mini')

    def prepare_host(self):
        # PDB file for the host, extract from the complex given as input.
        if not os.path.isfile(os.path.join(self.directory, self.guest + '.pdb')):
            subprocess.call('grep " ' + self.host + ' " ' + self.ipdb + '.pdb' + ' > ' + os.path.join(self.directory, self.guest + '.pdb'),
                            shell=True)
        if not os.path.isfile(os.path.join(self.directory, self.host + '.top')):
            self.generate_param_files(self.host)

    def minimize_complex(self):
        if not os.path.isfile(os.path.join(self.directory, 'topol-complex.top')):
            # create the topol-complex.top file
            f = open(os.path.join(self.static_dir, '00-modeling/complex.top'), 'r')
            filedata = f.read()
            f.close()

            newdata1 = filedata.replace('XXXXX', self.host)
            newdata2 = newdata1.replace('YYYYY', self.guest)

            f = open(os.path.join(self.directory, 'topol-complex.top'), 'w')
            f.write(newdata2)
            f.close()

        # First minimization step, steepest descent method.
        if not os.path.isfile(os.path.join(self.directory, 'mini1.trr')):
            self.minimize1(self.ipdb, 'topol-complex')

        # Second minimization step, conjugate gradient method.
        if not os.path.isfile(os.path.join(self.directory, 'mini2.trr')):
            self.minimize2('topol-complex')

        # Extract coordinates
        if not os.path.isfile(os.path.join(self.directory, 'complex_mini.pdb')):
            if not os.path.isfile(os.path.join(self.directory, 'TMP')):
                self.minimize2('topol-complex')
            self.extract_coordinates('complex_mini')

    def generate_param_files(self, who):
        subprocess.call('babel -ipdb ' + who + '.pdb -omol2 ' + who + '.mol2 --title ' + who,
                        shell=True)

        subprocess.call('cgenff ' + who + '.mol2 > ' + who + '.str',
                        shell=True)

        # Use cgenff_charmm2gmx (in the project) to create param for gmx
        mol_name = who
        mol2_name = who + '.mol2'
        rtp_name = who + '.str'
        GMXDATA = os.environ['GMXDATA']
        ffdir = os.path.join(GMXDATA, 'top/charmm36-jul2017.ff')
        cgenff_charmm2gmx.main(mol_name, mol2_name, rtp_name, ffdir)

        # Correct file names
        os.rename(who.lower() + '.top', who + '.top')
        os.rename(who.lower() + '.prm', who + '.prm')
        os.rename(who.lower() + '.itp', who + '.itp')
        os.rename(who.lower() + '_ini.pdb', who + '_ini.pdb')

    def minimize1(self, pdb, top):
        subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, '00-modeling/MINI1.mdp') + ' -c ' + pdb + '.pdb -p ' + top + '.top -o mini1.tpr -maxwarn 2',
                        shell=True)
        subprocess.call('gmx mdrun -v -deffnm mini1',
                        shell=True)

    def minimize2(self, top):
        subprocess.call('gmx_d grompp -f ' + os.path.join(self.static_dir, '00-modeling/MINI2.mdp') + ' -c mini1.gro -p ' + top + '.top -o mini2.tpr -maxwarn 2',
                        shell=True)
        subprocess.call('gmx mdrun -v -deffnm mini2 > TMP 2>&1',
                        shell=True)

    def extract_coordinates(self, pdb):
        # Find the number of the last step
        nstep = subprocess.check_output("grep 'Step' TMP | tail -1 | awk '{print $2}' | sed s/','//g",
                                        shell=True).decode("utf-8").rstrip('\n')
        # Extract the coordinate for the last step
        subprocess.call('echo "0" | gmx trjconv -f mini2.trr -s mini2.tpr -o ' + pdb + '.pdb -b ' + nstep,
                        shell=True)
        # clean
        os.remove(os.path.join(self.directory, 'TMP'))
