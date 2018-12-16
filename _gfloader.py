from neat import GreensTree

from _physiology import *

import numpy as np


def getGreensTree(morph_name, physiology_type='Pas', recompute=False, freqs=None):
    if freqs is None:
        freqs = np.array([0.])
    if morph_name[-4:] == '.swc':
        morph_name = morph_name[:-4]
    # load a greenstree
    greens_tree = GreensTree(file_n='morph/' + morph_name +'.swc')
    # set the physiological parameters
    eval('setPhysiology' + physiology_type + '(greens_tree)')
    # set the computational tree
    greens_tree.setCompTree(eps=1.)
    # set the impedance
    greens_tree.setImpedance(freqs)
    return greens_tree