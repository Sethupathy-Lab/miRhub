#!/usr/bin/env python

import os
import subprocess
import sys

class mirHubJob:

    def __init__(self):
        self.pro_net = 'miRhub/hsa_ppi_700_3col.txt'
        self.mjob = ['perl', 'miRhub/miRhub.pl']
        self.sjob = ['perl', 'miRhub/summary.pl']
        self.outname = []

    def add_DEG(self, x):
        if os.path.isfile(x):
            self.DEG = x
            file_base = os.path.basename(x)
            self.DEGbase = os.path.splitext(file_base)[0]
        else:
            print('ERROR: No gene list exists for provided name: {}'.format(x))
            sys.exit()
            
    def add_proj(self, x):
        if len(x) != len(x.strip()):
            print('ERROR: Remove spaces from project name')
            sys.exit()
        self.project = x

    def add_iter(self, x):
        self.iter = str(x)

    def add_cons(self, x):
        try:
            float(x)
        except ValueError:
            print('ERROR: cons argument contains non-numeric {}'.format(x))
            raise
        self.cons = x

    def add_spec_files(self, x):
        base = 'miRhub'
        if x == 'mouse':
            self.species = '10090'
            self.mirlist = '/'.join([base, 'all_mmu_mirs.txt'])
            self.scorecard = '/'.join([base, 'scorecard_MMU.txt'])
            os.path.isfile(self.scorecard)
        elif x == 'human':
            self.species = '9606'
            self.mirlist = '/'.join([base, 'all_hsa_mirs.txt'])
            self.scorecard = '/'.join([base, 'scorecard_HSA.txt'])
        elif x == 'rat':
            self.species = '10116'
            self.mirlist = '/'.join([base, 'all_rno_mirs.txt'])
            self.scorecard = '/'.join([base, 'scorecard_RNO.txt'])
        else:
            print('ERROR: {} not a valad species for miRhub'.format(x))

    def build_output_name(self):
        self.outname += [self.project,
                         self.DEGbase,
                         'cons{}'.format(self.cons),
                         'iter{}'.format(self.iter)]
        self.outname = '{}'.format('_'.join(self.outname))

    def build_mirhub_job(self):
        self.mjob += [self.mirlist,
                     self.DEG,
                     self.pro_net,
                     self.outname,
                     self.iter,
                     self.species,
                     self.cons]
        print(' '.join(self.mjob))
        os.system(' '.join(self.mjob))

    def build_summary_job(self):
        self.sjob += [self.outname,
                      self.DEG]
        os.system(' '.join(self.sjob))
        os.system('rm {}'.format(self.outname))
#        p = subprocess.Popen(self.sjob,
#                             stdout=subprocess.PIPE,
#                             stderr=subprocess.PIPE)
