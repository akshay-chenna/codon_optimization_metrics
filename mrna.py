from arnie.pfunc import pfunc
from arnie.free_energy import free_energy
from arnie.bpps import bpps
from arnie.mfe import mfe
import arnie.utils as utils

import rna_utils
import python_codon_tables as pct

import numpy as np
from decimal import Decimal

from draw_rna.ipynb_draw import draw_struct

import matplotlib.pyplot as plt
import seaborn as sns

from Bio.Seq import Seq
from Bio import SeqIO

import sys
sys.path.append('/home/akshay/apps/DegScore')

import os
os.environ['ETERNAFOLD_PATH'] = '/home/akshay/apps/miniconda3/envs/py311/bin/eternafold-bin/'
os.environ['ETERNAFOLD_PARAMETERS'] = '/home/akshay/apps/EternaFold/parameters/EternaFoldParams.v1'

from DegScore import DegScore

class rna_metrics:
    '''
    Get all mRNA metrics
    folding='eternafold' / 'vienna_2'
    seq = DNA sequence
    '''
    def __init__(self, fasta_file, folding='eternafold'):

        for record in SeqIO.parse(fasta_file, "fasta"):
            self.seq = record.seq

        self.fold = folding

    def arnie_free_energy(self):
        G = free_energy(self.seq, package=self.fold)
        return G

    def base_pair_matrix(self):
        bp_matrix = bpps(self.seq, package=self.fold)
        return bp_matrix

    def average_unpaired_probability(self):
        bp_matrix = bpps(self.seq, package=self.fold)
        p_unp_vec = 1 - np.sum(bp_matrix, axis=0)
        AUP = np.mean(p_unp_vec)
        return AUP
    
    def unpaired_probability(self):
        bp_matrix = bpps(self.seq, package=self.fold)
        p_unp_vec = 1 - np.sum(bp_matrix, axis=0)
        return p_unp_vec

    def mfe_structure(self):
        struct = mfe(self.seq,structure=self.fold)
        return struct

    def rna_2d_plots(self):
        struct = mfe(self.seq,structure=self.fold)
        draw_struct(self.seq, struct)

    def deg_score(self):
        mdl = DegScore(self.seq, package=self.fold)
        return mdl.degscore, mdl.est_half_life

    def CAI(self):
        dna = str(self.seq).replace("U","T")
        cc = rna_utils.CodonCollection(dna,pct.get_codons_table("h_sapiens_9606"))
        return cc.codon_adaptability_index()
