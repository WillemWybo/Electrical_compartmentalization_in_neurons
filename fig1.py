from neat import SOVTree, NET
import neat.netsim as netsim
from datarep import neuronmodel as neurm

from _matplotlibsettings import *
from _sovloader import *

import numpy as np

import sys
sys.setrecursionlimit(10000)


## paramteters #################################################################
morph_name = 'N19ttwt.CNG.swc'
locs = [{'node': 1, 'x': 0.5},
        {'node': 117, 'x': .9},
        {'node': 96, 'x': .3},
        {'node': 376, 'x':.5}]
v_eq = -75.
t_max = 300.
t0, t1, t2 = 10., 110., 210.
dt = .025
weight = 0.01
################################################################################


## deriving the NET ############################################################
# load the sov tree
sov_tree = getSOVTree(morph_name, physiology_type='Pas', recompute=False)
# derive the full NET
print '> Deriving the NET'
net_full = sov_tree.constructNET(dz=5., dx=10.,
                            add_lin_terms=False, improve_input_impedance=True)
# the points at which the impedance matrix is evaluated are stored under
# 'NET_eval', to prune the NET we need to find the indices of the locations in
# 'NET_eval' that are closes to our locations of interest
net_loc_inds = sov_tree.getNearestLocinds(locs, name='NET_eval')
# compute the impedance matrix of the full tree
z_xx = sov_tree.calcImpedanceMatrix(name='NET_eval')
# prune the full NET
net = net_full.getReducedTree(net_loc_inds, indexing='locs')
print net
################################################################################


## create and run the NEURON model #############################################
print '> Simulating NEURON model'
sim_tree = sov_tree.__copy__(new_tree=neurm.NeuronSimTree())
sim_tree.initModel(dt=dt, v_eq=v_eq, factor_lambda=10.)
# add synapses
sim_tree.addDoubleExpSynapse(locs[1], .2, 3., 0.)
sim_tree.addDoubleExpSynapse(locs[2], .2, 3., 0.)
sim_tree.addDoubleExpSynapse(locs[3], .2, 3., 0.)
# set the spike times
sim_tree.setSpikeTrain(0, weight, [t0])
sim_tree.setSpikeTrain(1, weight, [t1])
sim_tree.setSpikeTrain(2, weight, [t2])
# set recording locs
sim_tree.storeLocs(locs, name='rec locs')
# run the simulation
res_neuron = sim_tree.run(t_max)
################################################################################


## simulate the NET ############################################################
print '> Simulating the NET model'
cnet = netsim.NETSim(net, v_eq=v_eq)
# add synapses
cnet.addSynapse(1, "AMPA", g_max=weight)
cnet.addSynapse(2, "AMPA", g_max=weight)
cnet.addSynapse(3, "AMPA", g_max=weight)
# set the spike times
cnet.setSpikeTimes(0, np.array([t0])) # spike times for first synapse
cnet.setSpikeTimes(1, np.array([t1])) # spike times for second synapse
cnet.setSpikeTimes(2, np.array([t2])) # spike times for third synapse
# run the simulation
res_net = cnet.runSim(t_max, dt, rec_v_node=True)
################################################################################


## plot the figure #############################################################
print '> Plotting results'
fig = pl.figure('Fig1', figsize=(25,5.5))
set_title = True
spacings = [0.05,1.,0.4,1.,0.4,0.7,0.4,2.0,0.4,2.0,.1]
coords = np.cumsum(spacings)
coords /= coords[-1]

gs0_ = GridSpec(1, 1)
gs0_.update(top=.98, bottom=.85, left=coords[0], right=coords[1])
gs0 = GridSpec(1, 1)
gs0.update(top=.80, bottom=.15, left=coords[0], right=coords[1])

gs1_ = GridSpec(1, 1)
gs1_.update(top=.98, bottom=.85, left=coords[2], right=coords[3])
gs1 = GridSpec(1, 1)
gs1.update(top=.80, bottom=.15, left=coords[2], right=coords[3])

gs2_ = GridSpec(1, 1)
gs2_.update(top=.98, bottom=.85, left=coords[4], right=coords[5])
gs2 = GridSpec(1, 1)
gs2.update(top=.80, bottom=.15, left=coords[4], right=coords[5])

gs3_ = GridSpec(1, 1)
gs3_.update(top=.98, bottom=.85, left=coords[6], right=coords[7])
gs3 = GridSpec(4, 3)
gs3.update(top=.80, bottom=.15, left=coords[6], right=coords[7], wspace=.05, hspace=.05)

gs4_ = GridSpec(1, 1)
gs4_.update(top=.98, bottom=.85, left=coords[8], right=coords[9])
gs4 = GridSpec(5, 3)
gs4.update(top=.80, bottom=.15, left=coords[8], right=coords[9], wspace=.05, hspace=.05)

ax00 = pl.subplot(gs0_[0,0])
ax0 = pl.subplot(gs0[0,0])
ax00.axison = False
ax00.set_xlim((0.,0.))
ax00.set_ylim((0.,0.))
if set_title:
    ax00.text(0.,0., r'Morphology', fontsize=fontsize, horizontalalignment='center', verticalalignment='center')
sov_tree.plot2DMorphology(ax0,
                        plotargs={'color': 'DarkGrey', 'lw': 1.6},
                        marklocs=locs,
                        marklabels={0: r'soma', 1: '1', 2: '2', 3: '3'},
                        locargs=[{'marker': 'o', 'ms': 6, 'color': colours[ii]} for ii in range(len(locs))]
                       )

ax01 = pl.subplot(gs1_[0,0])
ax1 = pl.subplot(gs1[0,0])
ax01.axison = False
ax01.set_xlim((0.,0.))
ax01.set_ylim((0.,0.))
if set_title:
    ax01.text(0.,0., r'Impedance matrix', fontsize=fontsize, horizontalalignment='center', verticalalignment='center')
im = ax1.imshow(z_xx, origin='lower', interpolation='none', cmap='summer')
for i, ind in enumerate(net_loc_inds):
    ax1.axvline(ind, c=colours[i%len(colours)], lw=lwidth)
    ax1.axhline(ind, c=colours[i%len(colours)], lw=lwidth)
ax1.set_xlabel(r'$x$')
ax1.set_ylabel(r'$x^{\prime}$')
ax1.set_yticks([])
ax1.set_xticks([])
divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
cb = fig.colorbar(im, cax=cax)
cb.set_label(r'$Z(x,x^{\prime})$ (M$\Omega$)')

ax02 = pl.subplot(gs2_[0,0])
ax2 = pl.subplot(gs2[0,0])
ax02.axison = False
ax02.set_xlim((0.,0.))
ax02.set_ylim((0.,0.))
if set_title:
    ax02.text(0.,0., r'NET', fontsize=fontsize, horizontalalignment='center', verticalalignment='center')
syn_cs_d = {i: colours[i%len(colours)] for i in range(len(locs))}
syn_labels_d = {0: r'soma', 1: '1', 2: '2', 3: '3'}
# Zmax = net.plotDendrogram(ax2, reduced_tree, plotargs={'lw': lwidth}, textargs={'size': 'small'},
#                                     labelargs={'marker': 'o', 'ms': markersize, 'c': 'r'},
#                                     incolors=syn_cs_d, inlabels=syn_labels_d, nodelabels=None)

Zmax = net.plotDendrogram(ax2,
                            plotargs={'lw': lwidth},
                            textargs={'size': 'small'},
                            labelargs={'marker': 'o', 'ms': markersize, 'c': 'r'},
                            incolors=syn_cs_d,
                            inlabels=syn_labels_d,
                            nodelabels=None
                          )

ax04 = pl.subplot(gs3_[0,0])
ax04.axison = False
ax04.set_xlim((0.,0.))
ax04.set_ylim((0.,0.))
if set_title:
    ax04.text(0.,0., r'Membrane voltage', fontsize=fontsize, horizontalalignment='center', verticalalignment='center')
tplot = 100.
iplot = int(tplot/dt)
spiketimes1 = [t0, t1, t2]
spk_inds = [1,2,3]
for i, loc in enumerate(locs):
    labelstr = r'soma' if loc['node'] == 1 else r'$r_'+str(i)+'$'
    for j, spktm in enumerate(spiketimes1):
        i_start, i_stop = j*iplot, (j+1)*iplot
        # plot voltage
        ax_ = pl.subplot(gs3[i,j])
        # ax_.plot(res_net['t'][j*iplot:(j+1)*iplot]-j*tplot, res_net['v_loc'][i,j*iplot:(j+1)*iplot], ls='-',
        #                         c=colours[i%len(colours)], lw=lwidth, label=labelstr)
        # ax_.plot(res_neuron['t'][j*iplot:(j+1)*iplot]-j*tplot, res_neuron['v_m'][i,j*iplot:(j+1)*iplot], ls='--', c=colours[i%len(colours)], lw=lwidth*1.7)
        ax_.plot(res_net['t'][i_start:i_stop] - dt*i_start, res_net['v_loc'][i,i_start:i_stop], ls='-',
                                c=colours[i%len(colours)], lw=lwidth, label=labelstr)
        ax_.plot(res_neuron['t'][i_start:i_stop] - dt*i_start, res_neuron['v_m'][i,i_start:i_stop], ls='--', c=colours[i%len(colours)], lw=lwidth*1.7)
        # legend
        if j == len(spiketimes1)-1:
            leg = ax_.legend(loc=7, markerscale=lwidth, fontsize='small', borderpad=-1.0)
            # leg = ax_.legend(loc='upper right', markerscale=lwidth, fontsize='small')
            leg.draw_frame(False)
        # limits
        ylims = (-80.,-5.)
        ax_.set_ylim(ylims)
        ax_.set_xlim((0.,50.))
        # plot spike
        ax_.axvline(spktm-j*tplot, c=colours[j+1], lw=lwidth*2./3.)
        if i == 0:
            ax_.text(spktm+3.-j*tplot, ylims[0]+(ylims[1]-ylims[0])*0.9, r''+str(spk_inds[j]), fontsize='small')
        # labels
        if i == len(locs)-1 and j == 1:
            ax_.set_xlabel(r'$t$ (ms)')
        if j == 0 and i == 2:
            ax_.set_ylabel(r'$V_m$ (mV)')
        # ticks
        if i == len(locs)-1:
            ax_.xaxis.set_ticks([0.,20.,40.])
        else:
            ax_.xaxis.set_ticks([])
        if j == 0:
            ax_.yaxis.set_ticks([-80.,-40.])
        else:
            ax_.yaxis.set_ticks([])
        ax_.yaxis.set_ticks_position('left')
        ax_.xaxis.set_ticks_position('bottom')
        ax_.spines['top'].set_color('none')
        ax_.spines['right'].set_color('none')

ax03 = pl.subplot(gs4_[0,0])
ax03.axison = False
ax03.set_xlim((0.,0.))
ax03.set_ylim((0.,0.))
if set_title:
    ax03.text(0.,0., r'Node voltage', fontsize=fontsize, horizontalalignment='center', verticalalignment='center')
tplot = 100.
iplot = int(tplot/dt)
# i = 0
mapping = {(0,1,2,3): 4, (1,2): 3, (3,): 2, (2,): 1, (1,): 0}
for k, node in enumerate(net):
    if not (len(node.loc_inds) == 1 and 0 in node.loc_inds):
        node.loc_inds.sort()
        i = mapping[tuple(node.loc_inds)]
        # loop over spike times
        for j, spktm in enumerate(spiketimes1):
            # plot voltage
            ax_ = pl.subplot(gs4[i,j])
            if j == 2:
                ax_.plot(res_net['t'][j*iplot:(j+1)*iplot]-j*tplot, res_net['v_node'][k,j*iplot:(j+1)*iplot], c=cs[i],
                                                    lw=lwidth, label=r'$N=$'+''.join([str(ind) for ind in node.loc_inds]))
                leg = ax_.legend(loc=7, markerscale=lwidth, fontsize='small', borderpad=-2.0)
                leg.draw_frame(False)
            else:
                ax_.plot(res_net['t'][j*iplot:(j+1)*iplot]-j*tplot, res_net['v_node'][k,j*iplot:(j+1)*iplot], c=cs[i],
                                                    lw=lwidth)
            # limits
            ylims = (-5.,70.)
            ax_.set_ylim(ylims)
            ax_.set_xlim((0.,50.))
            # plot spike
            ax_.axvline(spktm-j*tplot, c=colours[j+1], lw=lwidth*2./3.)
            if i == 0:
                ax_.text(spktm+3.-j*tplot, ylims[0]+(ylims[1]-ylims[0])*0.9, r''+str(spk_inds[j]), fontsize='small')
            # labels
            if i == len(net)-2 and j == 1:
                ax_.set_xlabel(r'$t$ (ms)')
            if j == 0 and i == 2:
                ax_.set_ylabel(r'$\overline{V}_{N}$ (mV)')
            # ticks
            if i == len(net)-2:
                ax_.xaxis.set_ticks([0.,20.,40.])
            else:
                ax_.xaxis.set_ticks([])
            if j == 0:
                ax_.yaxis.set_ticks([0.,40.])
            else:
                ax_.yaxis.set_ticks([])
            ax_.yaxis.set_ticks_position('left')
            ax_.xaxis.set_ticks_position('bottom')
            ax_.spines['top'].set_color('none')
            ax_.spines['right'].set_color('none')

pl.show()
################################################################################
