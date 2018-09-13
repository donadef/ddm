# -*- coding: utf-8 -*-
import time
import configparser
import os
import numpy as np

from ddm.classes.modeling import Modeling
from ddm.classes.solvate import SolvateBound, SolvateUnbound
from ddm.classes.reference import PickReference
from ddm.classes.monitor import MonitorCVs
from ddm.classes.confine import ConfineBound, ConfineUnbound
from ddm.classes.vba import VbaBound,  VbaUnbound
from ddm.classes.alchemical import AlchemicalBound, AlchemicalUnbound
from ddm.classes.base import Guest, Host, Complex


class DDM:
    def __init__(self, config_file, pdb_complex):
        self.config = configparser.ConfigParser()
        self.config.read(config_file)

        self.pdb_complex = pdb_complex

        try:
            self.config['main']
        except configparser.NoSectionError:
            print("No [main] section in the config file")
            exit()

        try:
            self.host_name = self.config['main']['host']
        except KeyError:
            print("host parameter in [main] section of the config file is required.")
            exit()

        try:
            self.guest_name = self.config['main']['guest']
        except KeyError:
            print("guest parameter in [main] section of the config file is required.")
            exit()

        try:
            self.dest = self.config['main']['dest']
        except KeyError:
            print("dest parameter in [main] section of the config file is required.")
            exit()

        if not os.path.exists(self.dest):
            os.makedirs(self.dest)

        self.complex = Complex(self.pdb_complex, self.dest)
        self.guest = Guest(self.guest_name, self.complex, self.dest)
        self.host = Host(self.host_name, self.complex, self.dest)

    def perform_ddm(self):

        # First step : prepare the host, guest and complex.
        m = Modeling(self.config, self.guest, self.host, self.complex)
        m.run()

        # Solvate the bound and unbound states
        sb = SolvateBound(self.config, self.guest, self.host)
        sb.run()

        su = SolvateUnbound(self.config, self.guest)
        su.run()

        # Create reference and pick anchor points
        ref = PickReference(self.config)
        ref.run()

        # Collective variables
        mon = MonitorCVs(self.config)
        mon.run()

        # restrain LIGAND in BOUND state
        cf = ConfineBound(self.config, self.guest)
        dG_CONF_BND = cf.run()
        print('dG_CONF_BND', dG_CONF_BND)

        vba_bound = VbaBound(self.config, self.guest, mon.x0, mon.kappa_max, cf.krms_max)
        dG_VBA_BND = vba_bound.run()
        print('dG_VBA_BND', dG_VBA_BND)

        # decouple LIGAND in BOUND state
        alch_bound = AlchemicalBound(self.config, self.guest, mon.x0, mon.kappa_max, cf.krms_max)
        dG_ALCH_BND = alch_bound.run()
        print('dG_ALCH_BND', dG_ALCH_BND)

        # release LIGAND in UNBOUND state
        cf_unbound = ConfineUnbound(self.config, self.guest, cf.krms_max)
        dG_CONF_UB = cf_unbound.run()
        print('dG_CONF_UB', dG_CONF_UB)

        vba_unbound = VbaUnbound(self.config, mon.kappa_max)
        dG_VBA_UB, sym_corr = vba_unbound.run()
        print('dG_VBA_UB', dG_VBA_UB, 'sym_corr', sym_corr)

        # recouple LIGAND in UNBOUND state
        alch_unbound = AlchemicalUnbound(self.config, self.guest)
        dG_ALCH_UD = alch_unbound.run()
        print('dG_ALCH_UD', dG_ALCH_UD)

        all_contributions = dG_CONF_BND + dG_VBA_BND + dG_ALCH_BND + list(map(lambda x: -x, dG_CONF_UB)) + dG_VBA_UB + sym_corr + list(map(lambda x: -x, dG_ALCH_UD))
        total = np.sum(all_contributions)
        print('Total: ', total)
