from neat import SOVTree, NET
import neat.netsim as netsim
from datarep import neuronmodel as neurm

from _matplotlibsettings import *
from _sovloader import *
from _gfloader import *
from _spiketools import *

import copy

## parameters ##################################################################
morph_name = 'cell1_simplified.swc'
v_eq = -75. # mV
t_max = 1000. # ms
dt = .025 # ms
weights_exc = [0.00011, 0.00016, 0.00019]
# weights_exc = [0.0000, 0.,0.]
weights_inh = [0., 0.0023, 0.0046]
nmda_ratio = 3.
rate = 200. # Hz
rate_inh = 200. # Hz
# locations
soma_loc = {'node': 1, 'x': 0.5}
nmda_locs = [{'node': 306, 'x': 0.8}, {'node': 309, 'x': 0.8}] # synapses 3 and 4
# nmda_locs = [{'node': 246, 'x': 0.8}, {'node': 244, 'x': 0.8}] # synapses 1 and 2
shunt_loc = {'node': 303, 'x': 0.9}
locs = [soma_loc] + nmda_locs + [shunt_loc]
################################################################################


## deriving the NET ############################################################
# load the sov tree
sov_tree = getSOVTree(morph_name, physiology_type='ActiveSoma', recompute=False)
# derive the full NET
print '> Deriving the NET'
net_full = sov_tree.constructNET(dz=10., dx=10., eps=1e-3,
                            add_lin_terms=False, improve_input_impedance=True)
# the points at which the impedance matrix is evaluated are stored under
# 'NET_eval', to prune the NET we need to find the indices of the locations in
# 'NET_eval' that are closes to our locations of interest
net_loc_inds = sov_tree.getNearestLocinds(locs, name='NET_eval')
locs = [sov_tree.getLocs('NET_eval')[ind] for ind in net_loc_inds]
# prune the full NET
net = net_full.getReducedTree(net_loc_inds, indexing='locs')
################################################################################


################################################################################
greens_tree = getGreensTree(morph_name, physiology_type='ActiveSoma')
z_mat = greens_tree.calcImpedanceMatrix(locs).real[0]

z_bar_0123 = np.mean(z_mat[0,1:])
net[0].z_kernel.c *= z_bar_0123 / net[0].z_kernel.k_bar

z_bar_123 = np.mean([z_mat[1,3], z_mat[2,3], z_mat[1,2]]) - z_bar_0123
net[1].z_kernel.c *= z_bar_123 / net[1].z_kernel.k_bar

z_bar_12 = np.mean([z_mat[1,3], z_mat[2,3], z_mat[1,2]]) - z_bar_123 - z_bar_0123
net[2].z_kernel.c *= z_bar_12 / net[2].z_kernel.k_bar

z_bar_1 = z_mat[1,1] - z_bar_12 - z_bar_123 - z_bar_0123
net[3].z_kernel.c *= z_bar_1 / net[3].z_kernel.k_bar

z_bar_2 = z_mat[2,2] - z_bar_12 - z_bar_123 - z_bar_0123
net[4].z_kernel.c *= z_bar_1 / net[4].z_kernel.k_bar

z_bar_3 = z_mat[3,3] - z_bar_123 - z_bar_0123
net[5].z_kernel.c *= z_bar_3 / net[5].z_kernel.k_bar
################################################################################


reslist_net, reslist_neuron = [], []
cv_list, iz_list = [], []
for w_inh, w_exc in zip(weights_inh, weights_exc):
    ## spiketimes ##############################################################
    spktms_exc1 = poissonTrain(rate*1e-3, t_max)
    spktms_exc2 = poissonTrain(rate*1e-3, t_max)
    spktms_inh = np.arange(0., t_max, 1./(rate_inh*1e-3)).tolist()
    ############################################################################

    ## create and run the NEURON model #########################################
    print '> Simulating NEURON model'
    sim_tree = sov_tree.__copy__(new_tree=neurm.NeuronSimTree())
    sim_tree.treetype = 'computational'
    sim_tree.initModel(dt=dt, v_eq=v_eq, factor_lambda=10.)
    # add synapses
    sim_tree.addDoubleExpNMDASynapse(locs[1], .2, 3., .2, 43.,
                                     e_r=0., nmda_ratio=nmda_ratio)
    sim_tree.addDoubleExpNMDASynapse(locs[2], .2, 3., .2, 43.,
                                     e_r=0., nmda_ratio=nmda_ratio)
    sim_tree.addDoubleExpSynapse(locs[3], .2, 10., -80.)
    # set the spike times
    sim_tree.setSpikeTrain(0, w_exc, spktms_exc1) # spike times for nmda synapse 1
    sim_tree.setSpikeTrain(1, w_exc, spktms_exc2) # spike times for nmda synapse 2
    sim_tree.setSpikeTrain(2, w_inh, spktms_inh) # spike times for inhibitory synapse
    # set recording locs
    sim_tree.storeLocs(locs, name='rec locs')
    # run the simulation
    res_neuron = sim_tree.run(t_max)
    ############################################################################


    ## simulate the NET ########################################################
    print '> Simulating the NET model'
    cnet = netsim.NETSim(net, v_eq=v_eq)
    # add synapses
    cnet.addSynapse(1, 'AMPA+NMDA', g_max=w_exc, nmda_ratio=nmda_ratio)
    cnet.addSynapse(2, 'AMPA+NMDA', g_max=w_exc, nmda_ratio=nmda_ratio)
    cnet.addSynapse(3, 'GABA', g_max=w_inh)
    # set the spike times
    cnet.setSpikeTimes(0, spktms_exc1) # spike times for nmda synapse 1
    cnet.setSpikeTimes(1, spktms_exc2) # spike times for nmda synapse 2
    cnet.setSpikeTimes(2, spktms_inh) # spike times for inhibitory synapse
    # run the simulation
    res_net = cnet.runSim(t_max, dt, rec_v_node=True)
    ############################################################################

    ## compute correlations and Iz #############################################
    v_corr = np.corrcoef(res_net['v_loc'][1:3])[0,1]
    # compute average conductance GABA synapse
    tr, td = 0.2, 10.
    tp = (tr * td) / (td - tr) * np.log(td / tr)
    ff = 1. / (-np.exp(-tp / tr) + np.exp(-tp /td))
    s_gaba = ff * (td - tr)
    g_shunt = s_gaba * rate*1e-3 * w_inh
    # rescale net impedances according to Iz
    net_aux = copy.deepcopy(net)
    node_aux = net_aux.getLeafLocNode(1)
    z_bar = np.sum([node.z_kernel.k_bar for node in net_aux.pathToRoot(node_aux.parent_node)])
    for node in net_aux.pathToRoot(node_aux.parent_node):
        node.z_kernel.c /= (1. + z_bar * g_shunt)
    iz = net_aux.calcIZ([1,2])
    ############################################################################

    reslist_neuron.append(res_neuron); reslist_net.append(res_net)
    cv_list.append(v_corr); iz_list.append(iz)


## plot results ################################################################
pl.figure('v traces', figsize=(12,4))
axes = [pl.subplot(131), pl.subplot(132), pl.subplot(133)]
titels = [r'$w_{inh}$' + ' = %.5f uS'%(w_inh) for w_inh in weights_inh]
for res_neuron, res_net, titel, ax in zip(reslist_neuron, reslist_net, titels, axes):
    ax.set_title(titel)
    ax.plot(res_neuron['t'], res_neuron['v_m'][1], c='k')
    ax.plot(res_neuron['t'], res_neuron['v_m'][2], c='k')
    ax.plot(res_net['t'], res_net['v_loc'][1], c='r', ls='--', lw=lwidth)
    ax.plot(res_net['t'], res_net['v_loc'][2], c='b', ls='--', lw=lwidth)
    ax.set_xlabel(r'$t$ (ms)')
    ax.set_ylabel(r'$V_m$ (mV)')
    ax.set_ylim((-90.,0.))

pl.figure('Iz + vcorr', figsize=(8,4))
bar_width = 0.25
xx  = np.linspace(0.,1.,len(weights_inh))

ax = pl.subplot(121)
ax.bar(xx, cv_list, bar_width, color='b', alpha=.5, align='center')#, log=1)
ax.set_xlim((xx[0]-bar_width, xx[-1]+bar_width))
ax.set_xticks(xx)
ax.set_xticklabels(weights_inh)
ax.set_xlabel(r'$w_{inh}$ uS')
ax.set_ylabel(r'$C_{V_m}$')

ax = pl.subplot(122)
ax.bar(xx, iz_list, bar_width, color='b', alpha=.5, align='center')#, log=1)
ax.set_xlim((xx[0]-bar_width, xx[-1]+bar_width))
ax.set_xticks(xx)
ax.set_xticklabels(weights_inh)
ax.set_xlabel(r'$w_{inh}$ uS')
ax.set_ylabel(r'$I_Z$')

pl.tight_layout()
pl.show()
################################################################################

