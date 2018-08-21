# -*- coding: utf-8 -*-
import os
import shutil
import subprocess

class SolvateBound:
    def __init__(self, host, guest, complex, ipdb, dest, static_dir):
        self.host = host
        self.guest = guest
        self.complex = complex
        self.ipdb = ipdb
        self.dest = dest
        self.static_dir = static_dir
        self.prev_store = os.path.join(self.dest, '00-modeling/STORE')

        self.directory = os.path.join(self.dest, '01-solvate-bound')

    def run(self):
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

        os.chdir(self.directory)

        # Prepare simulation box
        if not os.path.isfile(os.path.join(self.directory, 'complex-cubic-box.pdb')):
            subprocess.call('echo "0" | gmx editconf -f ' + os.path.join(self.prev_store, 'complex_mini.pdb') + ' -o complex-cubic-box.pdb -bt cubic -d 1.2 -c -princ',
                            shell=True)

        # solvate
        if not os.path.isfile(os.path.join(self.directory, 'complex-cubic-box-solv.pdb')) or not os.path.isfile(os.path.join(self.directory, 'topol-complex-solv.top')) :
            subprocess.call('echo "0" | gmx solvate -cp complex-cubic-box.pdb -cs spc216.gro -o complex-cubic-box-solv.pdb > TMP 2>&1',
                            shell=True)

        nwat = subprocess.check_output("grep 'Number of SOL' TMP | awk '{print $5}'",
                               shell=True).decode("utf-8").rstrip('\n')

        # Generate restraints for solute
        if not os.path.isfile(os.path.join(self.directory, self.host + '-posre,itp')):
            subprocess.call('echo "0" | gmx genrestr -f ' + os.path.join(self.prev_store, self.host + '_ini.pdb') + ' -o ' + self.host + '-posre.itp',
                            shell=True)
        if not os.path.isfile(os.path.join(self.directory, self.guest + '-posre,itp')):
            subprocess.call('echo "0" | gmx genrestr -f ' + os.path.join(self.prev_store, self.host + '_ini.pdb') + ' -o ' + self.guest + '-posre.itp',
                            shell=True)

        # Create the topol file
        if not os.path.isfile(os.path.join(self.directory, 'topol-complex-solv.top')):
            f = open(os.path.join(self.static_dir, '01-solvate_bound/complex-solv.top'), 'r')
            filedata = f.read()
            f.close()

            newdata1 = filedata.replace('XXXXX', self.host)
            newdata2 = newdata1.replace('YYYYY', self.guest)
            newdata3 = newdata2.replace('ZZZZZ', nwat)

            f = open(os.path.join(self.directory, 'topol-complex.top'), 'w')
            f.write(newdata3)
            f.close()

        # Ionize
        if not os.path.isfile(os.path.join(self.directory, 'ions.tpr')):
            subprocess.call('gmx grompp -f MINI.mdp -c complex-cubic-box-solv.pdb -p topol-complex-solv.top -o ions.tpr -maxwarn 2 > TMP 2>&1',
                            shell=True)
        charge = subprocess.check_output("grep 'System has non-zero total charge' TMP | awk '{print $6}'",
                                         shell=True).decode("utf-8").rstrip('\n')

        if charge != '':  # TODO: charged complex
            if float(charge) > 0.0:
                pass
            else:
                pass
        else:
            shutil.copyfile('complex-cubic-box-solv.pdb', 'complex-cubic-box-solv-ions.pdb')

        # Run dynamics
        # MINI
        if not os.path.isfile(os.path.join(self.directory, 'mini.trr')):
            subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, '01-solvate_bound/MINI.mdp') + ' -c complex-cubic-box-solv-ions.pdb -p topol-complex-solv.top -o mini.tpr -maxwarn 2',
                            shell=True)
            subprocess.call('gmx mdrun -v -deffnm mini',
                            shell=True)
        # HEAT
        if not os.path.isfile(os.path.join(self.directory, 'heat.trr')):
            subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, '01-solvate_bound/HEAT.mdp') + ' -c mini.gro -p topol-complex-solv.top -o heat.tpr -maxwarn 2',
                            shell=True)
            subprocess.call('gmx mdrun -v -deffnm heat',
                            shell=True)

        # EQLB1
        if not os.path.isfile(os.path.join(self.directory, 'eqlb1.trr')):
            subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, '01-solvate_bound/EQLB1.mdp') + ' -c heat.gro -p topol-complex-solv.top -o eqlb1.tpr -maxwarn 2',
                            shell=True)
            subprocess.call('gmx mdrun -v -deffnm eqlb1',
                            shell=True)

        # EQLB2
        if not os.path.isfile(os.path.join(self.directory, 'eqlb2.trr')):
            subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, '01-solvate_bound/EQLB2.mdp') + ' -c eqlb1.gro -p topol-complex-solv.top -o eqlb2.tpr -maxwarn 2',
                            shell=True)
            subprocess.call('gmx mdrun -v -deffnm eqlb2',
                            shell=True)

        # PROD
        if not os.path.isfile(os.path.join(self.directory, 'prod.trr')):
            subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, '01-solvate_bound/PROD.mdp') + ' -c eqlb2.gro -p topol-complex-solv.top -o prod.tpr -maxwarn 2',
                            shell=True)
            subprocess.call('gmx mdrun -v -deffnm prod',
                            shell=True)


