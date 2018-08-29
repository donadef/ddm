# -*- coding: utf-8 -*-
import configparser
import os

from classes.modeling import Modeling
from classes.solvate import SolvateBound, SolvateUnbound
from classes.reference import PickReference
from classes.monitor import MonitorCVs
from classes.confine import ConfineBound, ConfineUnbound
from classes.vba import VbaBound,  VbaUnbound
from classes.alchemical import AlchemicalBound, AlchemicalUnbound
from classes.base import Guest, Host, Complex


class DDM:
    def __init__(self, config_file, pdb_complex):
        self.config = configparser.ConfigParser()
        self.config.read(config_file)

        self.pdb_complex = pdb_complex

        self.host_name = self.config['main']['host']
        self.guest_name = self.config['main']['guest']
        self.dest = self.config['main']['dest']
        if not os.path.exists(self.dest):
            os.makedirs(self.dest)

        self.complex = Complex(self.pdb_complex, self.dest)
        self.guest = Guest(self.guest_name, self.complex, self.dest)
        self.host = Host(self.host_name, self.complex, self.dest)

    def perform_ddm(self):
        # First step : prepare the host, guest and complex.

        m = Modeling(self.config, self.guest, self.host, self.complex)
        m.run()

        #
        sb = SolvateBound(self.config, self.guest, self.host)
        sb.run()

        su = SolvateUnbound(self.config, self.guest)
        su.run()

        ref = PickReference(self.config)
        ref.run()

        mon = MonitorCVs(self.config)
        mon.run()
        print(mon.x0)
        print(mon.kappa)
        print(mon.kappa_max)

        cf = ConfineBound(self.config, self.guest)
        dG_CONF_BND = cf.run()
        print(dG_CONF_BND)
        print(cf.krms)
        print(cf.krms_max)
        print(cf.flucts)

        vba_bound = VbaBound(self.config, self.guest, mon.x0, mon.kappa_max, cf.krms_max)
        dG_VBA_BND = vba_bound.run()
        print(dG_VBA_BND)

        vba_unbound = VbaUnbound(self.config, mon.kappa)
        dG_VBA_UB, sym_corr = vba_unbound.run()
        print(dG_VBA_UB, sym_corr)

        # alch_bound = AlchemicalBound(self.config, self.guest, mon.x0, mon.kappa_max, cf.krms_max)
        # dG_ALCH_BND = alch_bound.run()
        # print(dG_ALCH_BND)
        #
        # alch_unbound = AlchemicalUnbound(self.config, self.guest)
        # dG_ALCH_UD = alch_unbound.run()
        # print(dG_ALCH_UD)

        cf_unbound = ConfineUnbound(self.config, self.guest, cf.krms_max)
        cf_unbound.run()
