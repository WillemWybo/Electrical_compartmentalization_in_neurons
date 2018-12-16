import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as pl
import matplotlib.animation as manimation
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc, rcParams
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.offsetbox import AnchoredText

# colours = ['DeepPink', 'Purple', 'MediumSlateBlue', 'Blue', 'Teal',
#                 'ForestGreen',  'DarkOliveGreen', 'DarkGoldenRod',
#                 'DarkOrange', 'Coral', 'Red', 'Sienna', 'Black', 'DarkGrey']
colours = list(pl.rcParams['axes.prop_cycle'].by_key()['color'])


# matplotlib settings
legendsize = 10
labelsize = 15
ticksize = 15
lwidth = 1.5
markersize = 6.
fontsize = 16
lettersize = 20.
#~ font = {'family' : 'serif',
        #~ 'weight' : 'normal',
        #~ 'size'   : fontsize}
        #'sans-serif':'Helvetica'}
#'family':'serif','serif':['Palatino']}
#~ rc('font', **font)
rc('font',**{'family':'serif','serif':['Palatino'], 'size': 15.0})
rc('mathtext',**{'fontset': 'stixsans'})
# rc('text', usetex=True)
# rcParams['text.latex.preamble'].append(r"\usepackage{amsmath}\usepackage{xfrac}")
rc('legend',**{'fontsize': 'medium'})
rc('xtick',**{'labelsize': 'small'})
rc('ytick',**{'labelsize': 'small'})
rc('axes',**{'labelsize': 'large', 'labelweight': 'normal'})


cs = ['r', 'b', 'g', 'c', 'y']
mfs = ['D', 'o', 'v', '^', 's', 'p']
mls = ['+', '*', 'x', '1', '2']
lss = ['-', '--', '-.', ':']
cmap = pl.get_cmap('jet')

def myAx(ax):
    # customize the ax
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    return ax

def myLegend(ax, add_frame=True, **kwarg):
    leg = ax.legend(**kwarg)
    if add_frame:
        frame = leg.get_frame()
        frame.set_color('white')
        frame.set_alpha(0.8)
    return leg

def myColorbar(ax, im, **kwargs):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "5%", pad="3%")
    return pl.colorbar(im, cax=cax, **kwargs)


