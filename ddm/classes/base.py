# -*- coding: utf-8 -*-
import os
import shutil
import glob
import configparser


class DDMClass:
    def __init__(self, config, complex):
        self.config = config

        self.host = self.config['main']['host']
        self.guest = self.config['main']['guest']
        self.dest = self.config['main']['dest']

        self.complex = complex
        self.ipdb = os.path.basename(self.complex).rstrip('.pdb')

        self.static_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'static')
        self.awk_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'awk')

        self.directory = self.dest

    def store_files(self, files_to_store):
        store_dir = os.path.join(self.directory, 'STORE')
        if not os.path.exists(store_dir):
            os.makedirs(store_dir)

        for file in files_to_store:
            try:
                shutil.copy(file, store_dir)
            except FileNotFoundError:
                print('File ' + file + ' has not been created.')
                exit()

    def check_step(self, file):
        if not os.path.isfile(file):
            print("""
            ERROR : The process ended cause the file " + file + "has not been created as it should.
                    Please check the output above to find the error.""")
            exit()

    def clean_tmp(self):
        os.remove(os.path.join(self.directory, 'TMP'))

    def clean_md_files(self):
        for f in glob.glob("*.gro"):
            os.remove(f)
        for f in glob.glob("*.tpr"):
            os.remove(f)
        for f in glob.glob("*.trr"):
            os.remove(f)
        for f in glob.glob("*.edr"):
            os.remove(f)
        for f in glob.glob("*.log"):
            os.remove(f)
        for f in glob.glob("*.xtc"):
            os.remove(f)
        for f in glob.glob("*.mdp"):
            os.remove(f)
        for f in glob.glob("*#"):
            os.remove(f)


