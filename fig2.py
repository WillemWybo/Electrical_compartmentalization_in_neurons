from neat import SOVTree, NET
import neat.netsim as netsim
from datarep import neuronmodel as neurm

from _matplotlibsettings import *
from _sovloader import *
from _spiketools import *

import numpy as np
from neuronpy.util import spiketrain

import sys
sys.setrecursionlimit(10000)

## parameters ##################################################################
# morph_name = 'N19ttwt.CNG.swc'
# morph_name = 'stellate.swc'
morph_name = 'cell1_simplified.swc'
v_eq = -75. # mV
t_max = 10000. # ms
dt = .025 # ms
weight = 0.004 # uS
weight_inh = .5*weight # uS
Nsyn = 100
nmda_ratio = 3.
rate = 1. # Hz
rate_inh = rate / 2. # Hz
################################################################################


## spiketrains #################################################################
spiketimes = [poissonTrain(rate*1e-3, t_max) for _ in range(Nsyn+1)] + \
             [poissonTrain(rate_inh*1e-3, t_max) for _ in range(Nsyn+1)]
################################################################################


## deriving the NET ############################################################
# load the sov tree
sov_tree = getSOVTree(morph_name, physiology_type='ActiveSoma', recompute=False)
# derive the full NET
print '> Deriving the NET'
net_full, lin_terms_full = sov_tree.constructNET(dz=10., dx=10., eps=1e-3,
                            add_lin_terms=True, improve_input_impedance=True)
# define the input location indices
net_loc_inds = np.random.choice(np.arange(1,len(net_full[0].loc_inds)), Nsyn, replace=False)
net_loc_inds = np.concatenate(([0], net_loc_inds))
# prune the full NET
net = net_full.getReducedTree(net_loc_inds, indexing='locs')
# retain only necessary linear terms
lin_terms = {ii: lin_terms_full[ind] for ii,ind in enumerate(net_loc_inds) if ind != 0}
# rescale the weights if input impedance is very low
z_xx = net.calcImpMat()
if z_xx[0,0] < 80.:
    weight *= 2.
# compute somatic ion channel conductances
currents = {}
soma_surface = 4. * np.pi * (sov_tree[1].R*1e-4)**2 # cm^2
for channel_name, (g_bar, e_rev) in sov_tree[1].currents.iteritems():
    if channel_name != 'L':
        currents[channel_name] = (g_bar * soma_surface, e_rev)
################################################################################


## create and run the NEURON model #############################################
print '> Simulating NEURON model'
sim_tree = sov_tree.__copy__(new_tree=neurm.NeuronSimTree())
sim_tree.treetype = 'computational'
sim_tree.initModel(dt=dt, v_eq=v_eq, factor_lambda=10.)
# add synapses
locs = sov_tree.getLocs('NET_eval')
for loc_ind in net_loc_inds:
    sim_tree.addDoubleExpNMDASynapse(locs[loc_ind], .2, 3., .2, 43.,
                                     e_r=0., nmda_ratio=nmda_ratio)
for loc_ind in net_loc_inds:
    sim_tree.addDoubleExpSynapse(locs[loc_ind], .2, 10., -80.)
# set the spike times
for ii, loc_ind in enumerate(net_loc_inds):
    sim_tree.setSpikeTrain(ii, weight, spiketimes[ii])
n_exc = len(net_loc_inds)
for ii, loc_ind in enumerate(net_loc_inds):
    sim_tree.setSpikeTrain(ii+n_exc, weight_inh, spiketimes[ii+n_exc])
# set recording locs
rec_locs = [locs[ind] for ind in net_loc_inds]
sim_tree.storeLocs(rec_locs, name='rec locs')
# run the simulation
res_neuron = sim_tree.run(t_max, downsample=8, pprint=True)
################################################################################


## defining the simulation parameters ##########################################
print '> Simulating the NET model'
cnet = netsim.NETSim(net, lin_terms=lin_terms, v_eq=v_eq)
# cnet = netsim.NETSim(net, v_eq=v_eq)
# add the ion channels at the soma
for channel_name, (g_bar, e_rev) in currents.iteritems():
    cnet.addChannel(channel_name, 0, g_bar, e_rev=e_rev)
# # add the synapses
for loc_ind in range(len(net_loc_inds)):
    cnet.addSynapse(loc_ind, 'AMPA+NMDA', g_max=weight, nmda_ratio=nmda_ratio)
for loc_ind in range(len(net_loc_inds)):
    cnet.addSynapse(loc_ind, 'GABA', g_max=weight_inh, nmda_ratio=nmda_ratio)
# set the spike times
for ii, spktm in enumerate(spiketimes):
    cnet.setSpikeTimes(ii, spktm)
# run the simulation
res_net = cnet.runSim(t_max, dt, step_skip=8, pprint=True)
################################################################################


## compute membrane voltage correlations and spike coincidence #################
print '> analyzing results'
# compute correlation
v_corr_mat = np.corrcoef(res_net['v_loc'][1:])
i_z_mat = net.calcIZMatrix()[1:,1:]
v_corr_arr = np.reshape(v_corr_mat, v_corr_mat.size)
i_z_arr = np.reshape(i_z_mat, i_z_mat.size)
# get the spike times
spks_neuron = spikesFromVm(res_neuron['t'], res_neuron['v_m'][0])
spks_net = spikesFromVm(res_net['t'], res_net['v_loc'][0])
# get the coincidencefactor
if len(spks_neuron) > 0 and len(spks_net) > 0:
    cf = spiketrain.coincidence_factor(spks_neuron, spks_net, window=4.)
else:
    cf = np.nan
print '\n------ coincidence factor NET-NEURON spiketrains ------'
print 'cf =', cf
print '-------------------------------------------------------\n'
################################################################################


## plot results ################################################################
pl.figure('voltage')
ax = pl.gca()
ax.set_title('spike coincidence = %.2f'%cf)
ax.plot(res_neuron['t'], res_neuron['v_m'][0], 'r', label=r'NEURON')
ax.plot(res_net['t'], res_net['v_loc'][0], 'b--', label=r'NET')
ax.set_ylim((-90.,10.))
ax.set_xlabel(r'$t$ (ms)')
ax.set_ylabel(r'$V_m$ (mV)')
ax.legend(loc=0)

# to obtain a correlation plot as in the paper, simulating for at least
# 100 s is recommended
pl.figure('v correlations')
ax = pl.gca()
ax.scatter(i_z_arr, v_corr_arr)
ax.set_xlabel(r'$I_Z$')
ax.set_ylabel(r'$C_{V_m}$')

pl.show()
################################################################################
