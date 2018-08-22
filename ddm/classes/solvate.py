# -*- coding: utf-8 -*-
import os
import shutil
import subprocess

from .base import DDMClass


class Solvate(DDMClass):
    def __init__(self, config, complex):
        super(Solvate, self).__init__(config, complex)
        self.prev_store = os.path.join(self.dest, '00-modeling/STORE')

    def run(self):
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

        os.chdir(self.directory)

    def prepare_simulation_box(self, who):
        """
            Prepare the simulation box using Gromacs, editconf.

            -bt cubic
            -d 1.2
            -c -princ

            :param who: the 'complex' or the 'ligand' ?
            :type who: str
        """
        subprocess.call('echo "0" | gmx editconf -f ' + os.path.join(self.prev_store, who + '_mini.pdb') + ' -o ' + who + '-cubic-box.pdb -bt cubic -d 1.2 -c -princ',
                        shell=True)

    def solvate(self, who):
        """
            Solvate the solute using Gromacs, solvate.

            -cs spc216.gro

            :param who: the 'complex' or the 'ligand' ?
            :type who: str
            :return: The number of water added inside the box
            :rtype: str
        """
        subprocess.call('echo "0" | gmx solvate -cp ' + who + '-cubic-box.pdb -cs spc216.gro -o ' + who + '-cubic-box-solv.pdb > TMP 2>&1',
            shell=True)

        nwat = subprocess.check_output("grep 'Number of SOL' TMP | awk '{print $5}'",
                                       shell=True).decode("utf-8").rstrip('\n')
        self.clean_tmp()
        return nwat

    def generate_restraints(self, who):
        """
            Generate the restraints for the solute using Gromacs, genrestr.

            :param who: name of the solute (self.guest or self.host)
            :type who: str
        """
        subprocess.call('echo "0" | gmx genrestr -f ' + os.path.join(self.prev_store, who + '_ini.pdb') + ' -o ' + who + '-posre.itp',
                        shell=True)

    def ionize(self, who):
        """
            Ionize the system (if necessary) using Gromacs, genion.

            -pname NA
            -nname CL

            :param who: the 'complex' or the 'ligand' ?
            :type who: str
        """
        subprocess.call('gmx grompp -f MINI.mdp -c complex-cubic-box-solv.pdb -p topol-complex-solv.top -o ions.tpr -maxwarn 2 > TMP 2>&1',
                        shell=True)

        charge = subprocess.check_output("grep 'System has non-zero total charge' TMP | awk '{print $6}'",
                                         shell=True).decode("utf-8").rstrip('\n')
        self.clean_tmp()
        if charge != '':  # TODO: charged complex
            if float(charge) > 0.0:
                pass
            else:
                pass
        else:
            shutil.copyfile(who + '-cubic-box-solv.pdb', who + '-cubic-box-solv-ions.pdb')

    def run_dynamics(self, who):
        """
            Run dynamics using Gromacs, mdrun. 5 steps :

            -Minimization
            -Heat
            -Equilibration 1
            -Equilibration 2
            -Production

            :param who: the 'complex' or the 'ligand' ?
            :type who: str
        """
        # MINI
        if not os.path.isfile(os.path.join(self.directory, 'mini.trr')):
            subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, 'MINI.mdp') + ' -c ' + who + '-cubic-box-solv-ions.pdb -p topol-' + who + '-solv.top -o mini.tpr -maxwarn 2',
                            shell=True)
            subprocess.call('gmx mdrun -v -deffnm mini',
                            shell=True)
        # HEAT
        if not os.path.isfile(os.path.join(self.directory, 'heat.trr')):
            subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, 'HEAT.mdp') + ' -c mini.gro -p topol-' + who + '-solv.top -o heat.tpr -maxwarn 2',
                            shell=True)
            subprocess.call('gmx mdrun -v -deffnm heat',
                            shell=True)

        # EQLB1
        if not os.path.isfile(os.path.join(self.directory, 'eqlb1.trr')):
            subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, 'EQLB1.mdp') + ' -c heat.gro -p topol-' + who + '-solv.top -o eqlb1.tpr -maxwarn 2',
                            shell=True)
            subprocess.call('gmx mdrun -v -deffnm eqlb1',
                            shell=True)

        # EQLB2
        if not os.path.isfile(os.path.join(self.directory, 'eqlb2.trr')):
            subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, 'EQLB2.mdp') + ' -c eqlb1.gro -p topol-' + who + '-solv.top -o eqlb2.tpr -maxwarn 2',
                            shell=True)
            subprocess.call('gmx mdrun -v -deffnm eqlb2',
                            shell=True)

        # PROD
        subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, 'PROD.mdp') + ' -c eqlb2.gro -p topol-' + who + '-solv.top -o prod.tpr -maxwarn 2',
                        shell=True)
        subprocess.call('gmx mdrun -v -deffnm prod',
                        shell=True)


class SolvateBound(Solvate):
    def __init__(self, config, complex):
        super(SolvateBound, self).__init__(config, complex)
        self.directory = os.path.join(self.dest, '01-solvate-bound')
        self.static_dir = os.path.join(self.static_dir, '01-solvate-bound')

    def run(self):
        super(SolvateBound, self).run()

        files_to_store = []

        # Prepare simulation box
        if not os.path.isfile(os.path.join(self.directory, 'complex-cubic-box.pdb')):
            self.prepare_simulation_box('complex')

        # Solvate
        nb_wat = ''
        if not os.path.isfile(os.path.join(self.directory, 'complex-cubic-box-solv.pdb')) or not os.path.isfile(os.path.join(self.directory, 'topol-complex-solv.top')):
            nb_wat = self.solvate('complex')

        # Generate restraints for solute
        if not os.path.isfile(os.path.join(self.directory, self.host + '-posre.itp')):
            self.generate_restraints(self.host)
        if not os.path.isfile(os.path.join(self.directory, self.guest + '-posre.itp')):
            self.generate_restraints(self.guest)

        # Create the topol file
        if not os.path.isfile(os.path.join(self.directory, 'topol-complex-solv.top')) and not os.path.isfile(os.path.join(self.directory, 'STORE/topol-complex-solv.top')):
            if not os.path.isfile(os.path.join(self.directory, 'STORE/prod.tpr')):
                f = open(os.path.join(self.static_dir, 'complex-solv.top'), 'r')
                filedata = f.read()
                f.close()

                newdata = filedata.replace('XXXXX', self.host)
                newdata = newdata.replace('YYYYY', self.guest)
                newdata = newdata.replace('ZZZZZ', nb_wat)

                f = open(os.path.join(self.directory, 'topol-complex-solv.top'), 'w')
                f.write(newdata)
                f.close()
                files_to_store.append('topol-complex-solv.top')
        elif not os.path.isfile(os.path.join(self.directory, 'topol-complex-solv.top')) and os.path.isfile(os.path.join(self.directory, 'STORE/topol-complex-solv.top')):
            shutil.copy(os.path.join(self.directory, 'STORE/topol-complex-solv.top'), self.directory)

        # Ionize
        if not os.path.isfile(os.path.join(self.directory, 'complex-cubic-box-solv-ions.pdb')) and not os.path.isfile(os.path.join(self.directory, 'STORE/complex-cubic-box-solv-ions.pdb')):
            self.ionize('complex')
            files_to_store.append('complex-cubic-box-solv-ions.pdb')

        # Run dynamics
        if not os.path.isfile(os.path.join(self.directory, 'prod.tpr')) and not os.path.isfile(os.path.join(self.directory, 'STORE/prod.tpr')):
            self.run_dynamics('complex')
            files_to_store.append('prod.tpr')
            files_to_store.append('prod.xtc')

        self.store_files(files_to_store)

        self.clean_md_files()


class SolvateUnbound(Solvate):
    def __init__(self, config, complex):
        super(SolvateUnbound, self).__init__(config, complex)
        self.directory = os.path.join(self.dest, '02-solvate-unbound')
        self.static_dir = os.path.join(self.static_dir, '02-solvate-unbound')

    def run(self):
        super(SolvateUnbound, self).run()

        files_to_store = []

        # Prepare simulation box
        if not os.path.isfile(os.path.join(self.directory, 'ligand-cubic-box.pdb')):
            self.prepare_simulation_box('ligand')

        # Solvate
        nb_wat = ''
        if not os.path.isfile(os.path.join(self.directory, 'ligand-cubic-box-solv.pdb')) or not os.path.isfile(os.path.join(self.directory, 'topol-ligand-solv.top')):
            nb_wat = self.solvate('ligand')

        # Generate restraints for solute
        if not os.path.isfile(os.path.join(self.directory, self.guest + '-posre.itp')):
            self.generate_restraints(self.guest)

        # Create the topol file
        if not os.path.isfile(os.path.join(self.directory, 'topol-ligand-solv.top')) and not os.path.isfile(os.path.join(self.directory, 'STORE/topol-ligand-solv.top')):
            if not os.path.isfile(os.path.join(self.directory, 'STORE/prod.tpr')):
                f = open(os.path.join(self.static_dir, 'ligand-solv.top'), 'r')
                filedata = f.read()
                f.close()

                newdata = filedata.replace('YYYYY', self.guest)
                newdata = newdata.replace('ZZZZZ', nb_wat)

                f = open(os.path.join(self.directory, 'topol-ligand-solv.top'), 'w')
                f.write(newdata)
                f.close()
                files_to_store.append('topol-ligand-solv.top')
        elif not os.path.isfile(os.path.join(self.directory, 'topol-ligand-solv.top')) and os.path.isfile(os.path.join(self.directory, 'STORE/topol-ligand-solv.top')):
            shutil.copy(os.path.join(self.directory, 'STORE/topol-ligand-solv.top'), self.directory)

        # Ionize
        if not os.path.isfile(os.path.join(self.directory, 'ligand-cubic-box-solv-ions.pdb')) and not os.path.isfile(os.path.join(self.directory, 'STORE/ligand-cubic-box-solv-ions.pdb')):
            self.ionize('ligand')
            files_to_store.append('ligand-cubic-box-solv-ions.pdb')

        # Run dynamics
        if not os.path.isfile(os.path.join(self.directory, 'prod.trr')) and not os.path.isfile(os.path.join(self.directory, 'STORE/prod.tpr')):
            self.run_dynamics('ligand')
            files_to_store.append('prod.tpr')
            files_to_store.append('prod.xtc')

        self.store_files(files_to_store)

        self.clean_md_files()
