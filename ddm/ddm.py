# -*- coding: utf-8 -*-
import configparser
import os

from classes.modeling import Modeling
from classes.solvate import SolvateBound, SolvateUnbound
from classes.reference import PickReference
from classes.monitor import MonitorCVs
from classes.confine import ConfineBound
from classes.base import Guest


class DDM:
    def __init__(self, config_file, complex_file):
        self.config = configparser.ConfigParser()
        self.config.read(config_file)

        self.complex = complex_file

        self.host_name = self.config['main']['host']
        self.guest_name = self.config['main']['guest']
        self.dest = self.config['main']['dest']

        self.dihe = [0, 0, 1, 0, 1, 1]

    def perform_ddm(self):
        # First step : prepare the host, guest and complex.

        guest = Guest(self.guest_name, self.complex, self.dest)

        m = Modeling(self.config, self.complex, guest)
        m.run()

        sb = SolvateBound(self.config, self.complex)
        sb.run()

        su = SolvateUnbound(self.config, self.complex)
        su.run()

        ref = PickReference(self.config, self.complex)
        ref.run()

        mon = MonitorCVs(self.config, self.complex)
        mon.run()
        print(mon.x0)
        print(mon.kappa)
        print(mon.krms)

        cf = ConfineBound(self.config, self.complex)
        dG_CONF_BND = cf.run()
        print(cf.kappa)
        print(cf.krms)

