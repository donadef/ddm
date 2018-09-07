# -*- coding: utf-8 -*-
import os
import shutil
import subprocess
from alchemlyb.parsing.gmx import extract_u_nk
import pandas as pd
from alchemlyb.estimators import MBAR
from scipy import constants as c

from .base import DDMClass, check_step, clean_md_files, compute_std, compute_fluct, compute_trapez, compute_work
from .vba import Bound


ll_list = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23']


class AlchemicalBound(Bound):
    def __init__(self, config, guest, x0, kappa_max, krms_max):
        super(AlchemicalBound, self).__init__(config, guest, x0, kappa_max, krms_max)
        self.directory = os.path.join(self.dest, '08-alchemical-bound')
        self.static_dir = os.path.join(self.static_dir, '08-alchemical-bound')
        self.prev_store = os.path.join(self.dest, '06-vba-bound/STORE')
        self.prev_store_6 = os.path.join(self.dest, '05-confine-bound/STORE')
        self.prev_store_solv = os.path.join(self.dest, '01-solvate-bound/STORE')

        self.dG = []

    def run(self):
        super(AlchemicalBound, self).run()

        if not os.path.exists('STORE'):
            os.makedirs('STORE')

        if not os.path.isfile('dhdl-files/production-lambda_23.xvg'):
            if not os.path.exists('dhdl-files'):
                os.makedirs('dhdl-files')

            # Copy the topol file here
            shutil.copy(os.path.join(self.prev_store_solv, 'topol-complex-solv.top'), self.directory)

            plumedfile_data = self.create_plumed_tmpl()

            f = open(os.path.join(self.static_dir, 'PRODUCTION.tmpl'), 'r')
            prodfile_data = f.read()
            f.close()
            prodfile_data = prodfile_data.replace('XXXXX', self.guest.name)

            for ll in ll_list:
                if not os.path.isfile('dhdl-files/production-lambda_' + ll + '.xvg'):

                    new_plumedfile_data = plumedfile_data.replace('KK1', str(self.kappa_max[0]))
                    new_plumedfile_data = new_plumedfile_data.replace('KK2', str(self.kappa_max[1]))
                    new_plumedfile_data = new_plumedfile_data.replace('KK3', str(self.kappa_max[2]))
                    new_plumedfile_data = new_plumedfile_data.replace('KK4', str(self.kappa_max[3]))
                    new_plumedfile_data = new_plumedfile_data.replace('KK5', str(self.kappa_max[4]))
                    new_plumedfile_data = new_plumedfile_data.replace('KK6', str(self.kappa_max[5]))

                    f = open('plumed.dat', 'w')
                    f.write(new_plumedfile_data)
                    f.close()

                    new_prodfile_data = prodfile_data.replace('YYYYY', str(int(ll)))

                    f = open('PRODUCTION.mdp', 'w')
                    f.write(new_prodfile_data)
                    f.close()

                    subprocess.call('gmx grompp -f PRODUCTION.mdp -c ' + os.path.join(self.prev_store_6, '6.gro') + ' -t ' + os.path.join(self.prev_store_6, '6.cpt') + ' -p topol-complex-solv.top -o production.tpr -maxwarn 2',
                                    shell=True)
                    subprocess.call('gmx mdrun -deffnm production -plumed plumed.dat -v',
                                    shell=True)
                    check_step('production.xvg')
                    shutil.move('production.xvg', 'dhdl-files/production-lambda_' + ll + '.xvg')

                    clean_md_files()
            os.remove('plumed.dat')

        xvg_list = ['dhdl-files/production-lambda_' + ll + '.xvg' for ll in ll_list]
        u_nk = pd.concat([extract_u_nk(xvg, T=300) for xvg in xvg_list])
        mbar = MBAR().fit(u_nk)

        self.dG.append(mbar.delta_f_.iloc[0, len(xvg_list)-1] * c.Boltzmann * c.N_A / (4.184 * 1000) * 300)

        return self.dG


class AlchemicalUnbound(DDMClass):
    def __init__(self, config, guest):
        super(AlchemicalUnbound, self).__init__(config)
        self.directory = os.path.join(self.dest, '09-alchemical-unbound')
        self.static_dir = os.path.join(self.static_dir, '09-alchemical-unbound')
        self.prev_store = os.path.join(self.dest, '02-solvate-unbound/STORE')

        self.guest = guest

        self.dG = []

    def run(self):
        super(AlchemicalUnbound, self).run()

        if not os.path.exists('STORE'):
            os.makedirs('STORE')

        if not os.path.isfile('dhdl-files/production-lambda_23.xvg'):
            if not os.path.exists('dhdl-files'):
                os.makedirs('dhdl-files')

            # Copy the topol file here
            shutil.copy(os.path.join(self.prev_store, 'topol-ligand-solv.top'), self.directory)

            f = open(os.path.join(self.static_dir, 'PRODUCTION.tmpl'), 'r')
            prodfile_data = f.read()
            f.close()
            prodfile_data = prodfile_data.replace('XXXXX', self.guest.name)

            for ll in ll_list:
                if not os.path.isfile('dhdl-files/production-lambda_' + ll + '.xvg'):
                    new_prodfile_data = prodfile_data.replace('YYYYY', str(int(ll)))

                    f = open('PRODUCTION.mdp', 'w')
                    f.write(new_prodfile_data)
                    f.close()

                    subprocess.call('gmx grompp -f PRODUCTION.mdp -c ' + os.path.join(self.prev_store, 'prod.gro') + ' -t ' + os.path.join(self.prev_store, 'prod.cpt') + ' -p topol-ligand-solv.top -o production.tpr -maxwarn 2',
                                    shell=True)
                    subprocess.call('gmx mdrun -deffnm production -v',
                                    shell=True)
                    check_step('production.xvg')
                    shutil.move('production.xvg', 'dhdl-files/production-lambda_' + ll + '.xvg')

                    clean_md_files()

        xvg_list = ['dhdl-files/production-lambda_' + ll + '.xvg' for ll in ll_list]
        u_nk = pd.concat([extract_u_nk(xvg, T=300) for xvg in xvg_list])
        mbar = MBAR().fit(u_nk)

        self.dG.append(mbar.delta_f_.iloc[0, len(xvg_list)-1] * c.Boltzmann * c.N_A / (4.184 * 1000) * 300)

        return self.dG
