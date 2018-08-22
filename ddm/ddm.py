# -*- coding: utf-8 -*-
import configparser
import os

from classes.modeling import Modeling
from classes.solvate import SolvateBound, SolvateUnbound


class DDM:
    def __init__(self, config_file, complex_file):
        self.config = configparser.ConfigParser()
        self.config.read(config_file)

        self.complex = complex_file

    def perform_ddm(self):
        # First step : prepare the host, guest and complex.
        m = Modeling(self.config, self.complex)
        m.run()

        sb = SolvateBound(self.config, self.complex)
        sb.run()

        su = SolvateUnbound(self.config, self.complex)
        su.run()



