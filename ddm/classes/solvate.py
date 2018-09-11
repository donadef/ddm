# -*- coding: utf-8 -*-
import os
import shutil
import subprocess

from .base import DDMClass, clean_md_files, clean_tmp, check_step, ORGANIZE


class Solvate(DDMClass):
    def __init__(self, config):
        super(Solvate, self).__init__(config)
        self.prev_store = os.path.join(self.dest, ORGANIZE['modeling'], 'STORE')

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
        clean_tmp()
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
        clean_tmp()
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
        if not os.path.isfile('mini.gro'):
            subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, 'MINI.mdp') + ' -c ' + who + '-cubic-box-solv-ions.pdb -p topol-' + who + '-solv.top -o mini.tpr -maxwarn 2',
                            shell=True)
            subprocess.call('gmx mdrun -v -deffnm mini',
                            shell=True)
        # HEAT
        if not os.path.isfile('heat.gro'):
            subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, 'HEAT.mdp') + ' -c mini.gro -p topol-' + who + '-solv.top -o heat.tpr -maxwarn 2',
                            shell=True)
            subprocess.call('gmx mdrun -v -deffnm heat',
                            shell=True)

        # EQLB1
        if not os.path.isfile('eqlb1.gro'):
            subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, 'EQLB1.mdp') + ' -c heat.gro -p topol-' + who + '-solv.top -o eqlb1.tpr -maxwarn 2',
                            shell=True)
            subprocess.call('gmx mdrun -v -deffnm eqlb1',
                            shell=True)

        # EQLB2
        if not os.path.isfile('eqlb2.gro'):
            subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, 'EQLB2.mdp') + ' -c eqlb1.gro -p topol-' + who + '-solv.top -o eqlb2.tpr -maxwarn 2',
                            shell=True)
            subprocess.call('gmx mdrun -v -deffnm eqlb2',
                            shell=True)

        # PROD
        if not os.path.isfile('prod.gro'):
            subprocess.call('gmx grompp -f ' + os.path.join(self.static_dir, 'PROD.mdp') + ' -c eqlb2.gro -p topol-' + who + '-solv.top -o prod.tpr -maxwarn 2',
                            shell=True)
            subprocess.call('gmx mdrun -v -deffnm prod',
                            shell=True)


class SolvateBound(Solvate):
    def __init__(self, config, guest, host):
        super(SolvateBound, self).__init__(config)
        self.directory = os.path.join(self.dest, ORGANIZE['solvate-bound'])
        self.static_dir = os.path.join(self.static_dir, 'solvate-bound')

        self.guest = guest
        self.host = host

    def run(self):
        super(SolvateBound, self).run()

        # Prepare simulation box
        if not os.path.isfile('complex-cubic-box.pdb'):
            self.prepare_simulation_box('complex')

        # Solvate
        nb_wat = ''
        if not os.path.isfile('complex-cubic-box-solv.pdb') or not os.path.isfile('topol-complex-solv.top'):
            nb_wat = self.solvate('complex')

        # Generate restraints for solute
        if not os.path.isfile(self.host.name + '-posre.itp'):
            self.generate_restraints(self.host.name)
        if not os.path.isfile(self.guest.name + '-posre.itp'):
            self.generate_restraints(self.guest.name)

        # Create the topol file
        if not os.path.isfile('topol-complex-solv.top') and not os.path.isfile('STORE/topol-complex-solv.top'):
            f = open(os.path.join(self.static_dir, 'complex-solv.top'), 'r')
            filedata = f.read()
            f.close()

            newdata = filedata.replace('XXXXX', self.host.name)
            newdata = newdata.replace('YYYYY', self.guest.name)
            newdata = newdata.replace('ZZZZZ', nb_wat)
            if self.ff_param:
                newdata = newdata.replace("charmm36-jul2017.ff", self.ff_param)

            f = open('topol-complex-solv.top', 'w')
            f.write(newdata)
            f.close()

            self.store_files(['topol-complex-solv.top'])

        elif not os.path.isfile('topol-complex-solv.top') and os.path.isfile('STORE/topol-complex-solv.top'):
            shutil.copy('STORE/topol-complex-solv.top', self.directory)

        # Ionize
        if not os.path.isfile('complex-cubic-box-solv-ions.pdb') and not os.path.isfile('STORE/complex-cubic-box-solv-ions.pdb'):
            self.ionize('complex')
            self.store_files(['complex-cubic-box-solv-ions.pdb'])

        # Run dynamics
        if not os.path.isfile('STORE/prod.gro') or not os.path.isfile('STORE/prod.tpr'):
            self.run_dynamics('complex')
            check_step('prod.tpr')
            check_step('prod.xtc')
            check_step('prod.cpt')
            check_step('prod.gro')
            self.store_files(['prod.tpr', 'prod.xtc', 'prod.cpt', 'prod.gro'])

        clean_md_files()


class SolvateUnbound(Solvate):
    def __init__(self, config, guest):
        super(SolvateUnbound, self).__init__(config)
        self.directory = os.path.join(self.dest, ORGANIZE['solvate-unbound'])
        self.static_dir = os.path.join(self.static_dir, 'solvate-unbound')

        self.guest = guest

    def run(self):
        super(SolvateUnbound, self).run()

        # Prepare simulation box
        if not os.path.isfile('ligand-cubic-box.pdb'):
            self.prepare_simulation_box('ligand')

        # Solvate
        nb_wat = ''
        if not os.path.isfile('ligand-cubic-box-solv.pdb') or not os.path.isfile('topol-ligand-solv.top'):
            nb_wat = self.solvate('ligand')

        # Generate restraints for solute
        if not os.path.isfile(self.guest.name + '-posre.itp'):
            self.generate_restraints(self.guest.name)

        # Create the topol file
        if not os.path.isfile('topol-ligand-solv.top') and not os.path.isfile('STORE/topol-ligand-solv.top'):
            if not os.path.isfile('STORE/prod.tpr'):
                f = open(os.path.join(self.static_dir, 'ligand-solv.top'), 'r')
                filedata = f.read()
                f.close()

                newdata = filedata.replace('YYYYY', self.guest.name)
                newdata = newdata.replace('ZZZZZ', nb_wat)
                if self.ff_param:
                    newdata = newdata.replace("charmm36-jul2017.ff", self.ff_param)

                f = open(os.path.join(self.directory, 'topol-ligand-solv.top'), 'w')
                f.write(newdata)
                f.close()
                self.store_files(['topol-ligand-solv.top'])
        elif not os.path.isfile('topol-ligand-solv.top') and os.path.isfile('STORE/topol-ligand-solv.top'):
            shutil.copy('STORE/topol-ligand-solv.top', self.directory)

        # Ionize
        if not os.path.isfile('ligand-cubic-box-solv-ions.pdb') and not os.path.isfile('STORE/ligand-cubic-box-solv-ions.pdb'):
            self.ionize('ligand')
            self.store_files(['ligand-cubic-box-solv-ions.pdb'])

        # Run dynamics
        if not os.path.isfile('STORE/prod.gro') or not os.path.isfile('STORE/prod.tpr'):
            self.run_dynamics('ligand')
            check_step('prod.tpr')
            check_step('prod.xtc')
            check_step('prod.cpt')
            check_step('prod.gro')
            self.store_files(['prod.tpr', 'prod.xtc', 'prod.cpt', 'prod.gro'])

        clean_md_files()
