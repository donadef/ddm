# -*- coding: utf-8 -*-
import os
import shutil
import subprocess

from .base import DDMClass, check_step, clean_md_files, compute_std, compute_kf, compute_kf_plus, compute_fluct, compute_trapez


class Vba(DDMClass):
    def __init__(self, config, complex):
        super(Vba, self).__init__(config, complex)

    def run(self):
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

        os.chdir(self.directory)


class VbaBound(Vba):
    def __init__(self, config, complex):
        super(VbaBound, self).__init__(config, complex)
        self.directory = os.path.join(self.dest, '06-vba-bound')
        self.static_dir = os.path.join(self.static_dir, '06-vba-bound')
        self.prev_store = os.path.join(self.dest, '05-confine-bound/STORE')

    def run(self):
        super(VbaBound, self).run()


