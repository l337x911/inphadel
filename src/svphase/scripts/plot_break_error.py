import numpy as na
import matplotlib as mpl
from  matplotlib import pyplot as plt, gridspec
from operator import itemgetter

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100, N=100):
  new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
     'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
      cmap(na.linspace(minval, maxval, n)),N)
  return new_cmap

def plot_heatmap(h_r=12):
  npz = na.load('/home/anand/Projects/assembly/src/break_error_acc.npz')
  break_error_ranges = [(1,5),(5,10),(10,100)]

  fig = plt.figure(figsize=(11,6),facecolor=None, edgecolor=None)
  
  gs = gridspec.GridSpec(2,3, height_ratios=[h_r,1])
  #gs.update(left=0.05, bottom=0.05, top=0.95, right=0.95, wspace=0.01, hspace=0.01)

  z_min, z_max = 0,1
  #cmap = truncate_colormap(plt.cm.YlGnBu,0.2,0.95,N=50)
  cmap = plt.cm.YlGnBu
   
  for idx, (r, a) in enumerate(zip(break_error_ranges, map(itemgetter(1), sorted(npz.items())))):
    ax = plt.subplot(gs[0,idx], aspect='equal')
    offsets = map('{0:d}'.format, na.arange(0,a.shape[0]*100,100))
    hm = ax.pcolor(a, cmap=cmap, vmin=z_min, vmax=z_max)

    # put the major ticks at the middle of each cell
    ax.set_xticks(na.arange(a.shape[0])+0.5, minor=False)
    ax.set_yticks(na.arange(a.shape[1])+0.5, minor=False)
    
    ax.set_xticklabels(offsets, minor=False)
    ax.set_yticklabels(offsets, minor=False)
    ax.set_xlabel('Deletions with size\n{0}kb-{1}kb'.format(*r))
    ax.set_ylabel('Offset from true break a')

    ax.invert_yaxis()
    ax.xaxis.tick_top()    
  
  cax =  plt.subplot(gs[1,:])
  cb = fig.colorbar(hm, cax=cax, orientation='horizontal', format='%.2f')  
  cax.set_xlabel('accuracy')
  fig.tight_layout()
  plt.show()

if __name__ == '__main__':
  plot_heatmap()
