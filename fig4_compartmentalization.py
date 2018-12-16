from neat import SOVTree, NET
import neat.netsim as netsim
from datarep import neuronmodel as neurm

from _matplotlibsettings import *
from _sovloader import *

import numpy as np

## parameters ##################################################################
# morph_name = 'N19ttwt.CNG.swc'
# morph_name = 'stellate.swc'
morph_name = 'cell1_simplified.swc'
Iz = 10.
################################################################################


## deriving the NET ############################################################
# load the sov tree
sov_tree = getSOVTree(morph_name, physiology_type='ActiveSoma', recompute=False)
# derive the full NET
print '> Deriving the NET'
net = sov_tree.constructNET(dz=20., dx=10., add_lin_terms=False)
# get the compartments
print '> Deriving compartmentalization'
cnode_inds = net.getCompartmentalization(Iz)
################################################################################


## plot the compartments #######################################################
print '> Plotting results'
n_comp = len(cnode_inds)
comp_colors = np.random.permutation(np.linspace(0., 1., n_comp))
# color dictionary for NET dendrogram
cs_comp_net = {ind: comp_colors[ii] for ii, c_inds in enumerate(cnode_inds) for ind in c_inds}
# color dictionary for morphology
cs_comp_morph = {}
# convert net nodes to morph nodes
for inds in cnode_inds:
    for ind in inds:
        loc_inds = net[ind].loc_inds
        node_inds = sov_tree.getNodeIndices('NET_eval')[loc_inds]
        ii = np.argmin(sov_tree.distancesToSoma('NET_eval')[loc_inds])
        cs_comp_morph.update({node.index: cs_comp_net[ind] for node in sov_tree.__iter__(node=sov_tree[node_inds[ii]])})
# plot the figure
pl.figure()
ax0 = pl.subplot(121)
sov_tree.plot2DMorphology(ax0, cs=cs_comp_morph, use_radius=False, cmap=cmap,
                          plotargs={'color': 'DarkGrey', 'lw': lwidth})
ax1 = pl.subplot(122)
net.plotDendrogram(ax1, cs_comp=cs_comp_net, cmap=cmap,
                    plotargs={'color': 'DarkGrey', 'lw': lwidth})
pl.show()
################################################################################

