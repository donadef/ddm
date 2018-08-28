# -*- coding: utf-8 -*-
import configparser
import os

from classes.modeling import Modeling
from classes.solvate import SolvateBound, SolvateUnbound
from classes.reference import PickReference
from classes.monitor import MonitorCVs
from classes.confine import ConfineBound
from classes.vba import VbaBound,  VbaUnbound
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

    def perform_ddm(self):
        # First step : prepare the host, guest and complex.

        complex = Complex(self.pdb_complex, self.dest)
        guest = Guest(self.guest_name, complex, self.dest)
        print(guest.beg)
        print(guest.end)

        host = Host(self.host_name, complex, self.dest)

        m = Modeling(self.config, guest, host, complex)
        m.run()
        #
        sb = SolvateBound(self.config, guest, host)
        sb.run()

        su = SolvateUnbound(self.config, guest)
        su.run()

        ref = PickReference(self.config)
        ref.run()

        mon = MonitorCVs(self.config)
        mon.run()
        print(mon.x0)
        print(mon.kappa)
        print(mon.kappa_max)

        cf = ConfineBound(self.config, guest)
        dG_CONF_BND = cf.run()
        print(cf.krms)
        print(cf.krms_max)
        print(cf.flucts)

        vba_bound = VbaBound(self.config, guest, mon.x0, mon.kappa_max, cf.krms_max)
        dG_VBA_BND = vba_bound.run()

        print(dG_VBA_BND)

        vba_unbound = VbaUnbound(self.config, mon.kappa_max)
        dG_VBA_UB = vba_unbound.run()


