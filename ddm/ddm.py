# -*- coding: utf-8 -*-
import configparser
import os
from classes.modeling import Modeling
from classes.solvate_bound import SolvateBound

class DDM:
    def __init__(self, config, complex):
        self.config = configparser.ConfigParser()
        self.config.read(config)

        self.host = self.config['main']['host']
        self.guest = self.config['main']['guest']
        self.dest = self.config['main']['dest']

        self.complex = complex
        self.ipdb = os.path.basename(self.complex).rstrip('.pdb')

        self.static_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'ddm/static')

    def perform_ddm(self):
        # First step : prepare the host, guest and complex.
        # m = Modeling(self.host, self.guest, self.complex, self.ipdb, self.dest, self.static_dir)
        # m.run()

        sb = SolvateBound(self.host, self.guest, self.complex, self.ipdb, self.dest, self.static_dir)
        sb.run()
