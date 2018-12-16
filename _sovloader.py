from neat import SOVTree

from _physiology import *

import dill
dill.settings['recurse'] = True


def getSOVTree(morph_name, physiology_type='Pas', recompute=False):
    print '> Computing SOV'
    if morph_name[-4:] == '.swc':
        morph_name = morph_name[:-4]
    file_name = morph_name + '_' + physiology_type + '.p'
    try:
        sov_tree = dill.load(open('sov_models/' + file_name, 'rb'))
        # ensure that the tree is recomputed if 'recompute' is true
        if recompute:
            raise IOError
    except (IOError, EOFError) as err:
        print '\'' + file_name + '\' does not exist, calculating cell...'
        sov_tree = SOVTree(file_n='morph/' + morph_name +'.swc')
        # set the physiological parameters
        eval('setPhysiology' + physiology_type + '(sov_tree)')
        # set the computational tree
        sov_tree.setCompTree(eps=1.)
        # compute SOV factorisation
        sov_tree.calcSOVEquations()
        # store the tree
        dill.dump(sov_tree, open('sov_models/' + file_name, 'wb'))

    return sov_tree

