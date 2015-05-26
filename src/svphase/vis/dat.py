'''
Created on Jan 29, 2014

@author: anand
'''
from matplotlib import pyplot as plt, gridspec
import matplotlib as mpl
import numpy as na
from svphase.parser import ReadsParserDat
from itertools import izip 

class TriadColors(object):
  def __init__(self, n=100):
    hexA,hexB,hexC = '#0011FF','#FF9C00','#40B209'
    hexR = '#FF0000'
    white = (1.0, 1.0, 1.0, 0.0)

    self.cmaps = [mpl.colors.LinearSegmentedColormap.from_list('linBlue', [white, hexA], n),
                  mpl.colors.LinearSegmentedColormap.from_list('linOrange', [white, hexB], n),
                  mpl.colors.LinearSegmentedColormap.from_list('linGreen', [white, hexC], n),
                  mpl.colors.LinearSegmentedColormap.from_list('linRed', [white, hexR], n)]
    self.cmaps = [truncate_cmap(cm, 0.2, 1, n=n) for cm in self.cmaps]
    n = self.cmaps[0].N
    #alphas = na.zeros(n+3)
    alphas = na.ones(n+3)
    #n_s = 1*n/2
    #alphas[:n_s] = na.linspace(.6,1,n_s)
    #alphas[n_s:n] = 1
    for hc, cm in izip([hexA, hexB, hexC, hexR],self.cmaps):
      cm._init()
      cm._lut[:,-1] = alphas
      cm.set_under(color='w',alpha=0)
      cm.set_over(color=hc, alpha=1)
      cm.set_bad(alpha=0)
def truncate_cmap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(na.linspace(minval, maxval, n)))
    return new_cmap
  
class VisualizeDat(object):
  def __init__(self, fpath):
    self.reads_dat = ReadsParserDat()
    self.reads_dat._index(fpath)
    self.current_fpath = fpath
  def _set_max_p(self, resolution):
    self.max_p = 1+na.max(self.reads_dat.index)/resolution
  def _load_matrices(self, contig, resolution):
    self._set_max_p(resolution)
    
    verted = na.zeros((self.max_p, self.max_p), dtype=na.uint16)
    fr = na.zeros((self.max_p,self.max_p), dtype=na.uint16)
    rf = na.zeros((self.max_p,self.max_p), dtype=na.uint16)

    for row in self.reads_dat.get_reads(self.current_fpath, contig):
      posa,posb = row[1]/resolution,row[3]/resolution
      strda,strdb = self.reads_dat.get_strandedness(row)
      if strda==strdb:
        verted[posb,posa] += 1
      elif strda=='+':
        fr[posb,posa] += 1
      else:
        rf[posb,posa] += 1
    return verted, fr, rf    

  def _load_matrices_avg(self, contig, resolution, window_size=100, n=100):
    self._set_max_p(resolution)
    
    verted = na.zeros((window_size, window_size), dtype=na.uint16)
    fr = na.zeros((window_size,window_size), dtype=na.uint16)
    rf = na.zeros((window_size,window_size), dtype=na.uint16)
    
    p = na.random.randint(0,self.max_p-window_size+1, size=n)
    for row in self.reads_dat.get_reads(self.current_fpath, contig):
      posa,posb = row[1]/resolution,row[3]/resolution
      strda,strdb = self.reads_dat.get_strandedness(row)
      offa,offb = posa-p, posb-p
      offa_m =  na.logical_and(offa>=0, offa<window_size)

      if not offa_m.any(): continue

      offb_m =  na.logical_and(offb>=0, offb<window_size)
      window_idx = na.nonzero(na.logical_and(offa_m,offb_m))[0]

      if len(window_idx)==0: continue
      # assumes non-overlapping windows.
      window_idx = window_idx[0]
      
      nposa,nposb = offa[window_idx],offb[window_idx]
      assert nposa>=0
      assert nposa<window_size
      assert nposb>=0
      assert nposb<window_size
     
      if strda==strdb:
        verted[nposb,nposa] += 1
      elif strda=='+':
        fr[nposb,nposa] += 1
        assert nposa<=nposb
      else:
        rf[nposb,nposa] += 1
    return verted, fr, rf 

  def visualize(self, contig, ax, resolution=100000, small_diag=None, norm_const=None, color_ax=None,vmax=500, avg_flag=False, single_color=False):
    self._set_max_p(resolution)

    if small_diag is None:
      verted, fr, rf = self._load_matrices(contig, resolution)
      vmin = 1
      norm_func = mpl.colors.LogNorm
    elif avg_flag:
      verted, fr, rf = map(lambda x:na.array(x,dtype=na.float)/small_diag[1], self._load_matrices_avg(contig, resolution,small_diag[0], small_diag[1]))
      vmin = 1./(small_diag[1])
      vmax = 4.0
      norm_func = mpl.colors.Normalize
    else:
      verted, fr, rf = self._load_matrices_avg(contig, resolution,small_diag[0], small_diag[1])
      vmin = 1
      vmax = 15
      norm_func = mpl.colors.Normalize
      #norm_func = mpl.colors.LogNorm

    #na.savez('s.npz', v=verted, fr=fr, rf=rf)
    #tmp_array = na.load('t.npz')
    #verted,fr,rf = [tmp_array[i] for i in ['v', 'fr', 'rf']]

    if not norm_const is None:
      vmax=20
      verted = na.array(verted,dtype=na.float)/norm_const
      fr = na.array(fr,dtype=na.float)/norm_const
      rf = na.array(rf,dtype=na.float)/norm_const
      vmin = 0.001

    if len(ax)==1:
      ax1 = ax[0]
      ax2 = ax[0]
      ax3 = ax[0]
      alpha = 0.5
    else:
      ax1,ax2,ax3 = ax
      alpha = 1

    bins =na.concatenate(([0],[norm_func(vmin=vmin,vmax=vmax).inverse(i) for i in na.linspace(0,1,20)]))
    print bins
    print na.max(verted), na.sum(verted), na.histogram(verted[verted>0],bins=bins)[0]
    print na.max(fr), na.sum(fr), na.histogram(fr[fr>0], bins=bins)[0]
    print na.max(rf), na.sum(rf), na.histogram(rf[rf>0], bins=bins)[0]
    print "Window Size", self.max_p, "Max Value", vmax
    triad = TriadColors()

    ax_cmaps = [triad.cmaps[2],
                  #cmap=truncate_cmap(plt.get_cmap('Greens'), 0.3,1.0,500),
                  triad.cmaps[1],
                  #cmap=truncate_cmap(plt.get_cmap('Blues'), 0.3, 1.0,500),
                  triad.cmaps[0],
                  #cmap=truncate_cmap(plt.get_cmap('Reds'), 0.3, 1.0,500),
                  triad.cmaps[3]]

    for ax_num, single_ax in enumerate(ax):
      single_ax.pcolormesh(na.zeros_like(verted), cmap=plt.get_cmap('binary'), edgecolor='None', lw=0)

    total = (verted+fr+rf).T
    labels_to_cmap_dict = dict(zip(['-/-,+/+', '-/+', '+/-', 'all'], range(4)))

    if single_color:
      for ax_num, single_ax in enumerate(ax):
        single_ax.pcolormesh(total.T,cmap=ax_cmaps[3],
                              norm=norm_func(vmin=vmin, vmax=vmax), edgecolor='None', lw=0)
    else:  
      ax3.pcolormesh(verted, cmap=ax_cmaps[0],
                    norm=norm_func(vmin=vmin, vmax=vmax), edgecolor='None', lw=0, alpha=alpha)
      ax2.pcolormesh(rf,cmap=ax_cmaps[1],
                    norm=norm_func(vmin=vmin, vmax=vmax), edgecolor='None', lw=0, alpha=alpha)
      ax1.pcolormesh(fr,cmap=ax_cmaps[2],
                    norm=norm_func(vmin=vmin, vmax=vmax), edgecolor='None', lw=0, alpha=alpha)

    for ax_num, single_ax in enumerate(ax):
      single_ax.pcolormesh(total,cmap=ax_cmaps[3],
                  norm=norm_func(vmin=vmin, vmax=vmax), edgecolor='None', lw=0)
      #single_ax.invert_yaxis()
      single_ax.set_aspect('equal')
      single_ax.axis('off')

    ax1.axis('on')
    ax1.set_frame_on(False)
    ax1.get_xaxis().tick_top()
    ax1.get_yaxis().tick_left()

    if not color_ax is None:
      vspace = 0.02
      cx = na.linspace(vmin,vmax*3,1024)
      cy = [0,1]
      cX, cY = na.meshgrid(cx,cy)
      
      if single_color:
        displays = ['all'] 
      else:
        displays = ['-/-,+/+', '-/+', '+/-', 'all'] 

      for i, label  in enumerate(displays):
        norm = norm_func(vmin=vmin, vmax=vmax)
        C = cx.reshape(1,-1)
        C = na.vstack((C,C))
        color_ax.pcolormesh(cX,cY+i+i*vspace,C, cmap=ax_cmaps[labels_to_cmap_dict[label]], norm=norm, edgecolors='None')
        color_ax.text(vmin,0.5+i+i*vspace,label, verticalalignment='center')
        #for rx in lin_v[::8]:
        #  color_ax.text(norm(rx),0.75+i+i*vspace, "%.2f"%rx, ha='center') 
      color_ax.get_xaxis().tick_bottom()
      color_ax.axes.get_yaxis().set_visible(False)
      color_ax.set_yticks([])
      if small_diag is None:
        color_ax.set_xscale('log')
      
      color_ax.set_xlim(vmin,vmax) 

      color_ax.set_xlabel('Number of Reads') 

  def visualize_fig(self, contig, resolution=100000, w_r=5, small_diag=None, separate_flag=False):

    if separate_flag:
      gs = gridspec.GridSpec(3,2, width_ratios=[1,1], height_ratios=[w_r,w_r,1])
    else:
      gs = gridspec.GridSpec(2,1, height_ratios=[w_r,1])
    gs.update(left=0.05, bottom=0.05, top=0.95, right=0.95, wspace=0.01, hspace=0.01)

    ax = plt.subplot(gs[0,0])
    color_ax = plt.subplot(gs[-1,:])

    if separate_flag:
      ax2 = plt.subplot(gs[0,-1], sharex=ax,sharey=ax)
      ax3 = plt.subplot(gs[-2,0], sharex=ax,sharey=ax)
      self.visualize(contig, [ax,ax2,ax3], resolution, small_diag=small_diag, color_ax=color_ax)
    else:
      self.visualize(contig, [ax,], resolution, small_diag=small_diag,  color_ax=color_ax)

    if small_diag is None:
      self._set_max_p(resolution)
      max_p = self.max_p
    else:
      max_p = small_diag[0]
    
    ax.set_xlim(0,max_p)
    ax.set_ylim(0,max_p)
    ax.invert_yaxis()

  def _format_label(self, ax, resolution, small_diag):
    if small_diag is None:
      # use MB
      r = 1000000/resolution
      label_func = lambda x,l: "{0:.1f} MB".format(x/r)
      label_size = 'small'
      label_size = 'medium'
    else:
      # use KB
      r = 1000/resolution
      label_func = lambda x,l: "{0:.1f} KB".format(x/r)
      #label_size = 'x-large'
      label_size = 'large'
    #print resolution, r, small_diag
    #print ax.get_xticks()
    #print ax.get_yticks()
    #ax.set_xticklabels( map(label_func, ax.get_xticks()), size=label_size)
    #ax.set_yticklabels( map(label_func, ax.get_yticks()), size=label_size)
    fmt = mpl.ticker.FuncFormatter(label_func)
    ax.xaxis.set_major_formatter(fmt)
    ax.yaxis.set_major_formatter(fmt)
    ax.tick_params(labelsize=label_size)

 
  def visualize_linear_fig(self, contig, resolution=100000, w_r=5, small_diag=None):
    fig = plt.figure(figsize=(14,6), facecolor='white', edgecolor='none', frameon=False)
    #fig = plt.figure(figsize=(14,6), facecolor='white', edgecolor='none', frameon=True)
    gs = gridspec.GridSpec(2,3, width_ratios=[1,1,1], height_ratios=[w_r,1])
    gs.update(left=0.15, bottom=0.10, top=0.95, right=0.95, wspace=0, hspace=0.01)

    ax = fig.add_subplot(gs[0,0])
    color_ax = fig.add_subplot(gs[-1,:])
    color_ax.set_frame_on(False)

    ax2 = plt.subplot(gs[0,1], sharex=ax,sharey=ax)
    ax3 = plt.subplot(gs[0,2], sharex=ax,sharey=ax)
    self.visualize(contig, [ax,ax2,ax3], resolution, small_diag=small_diag, color_ax=color_ax)

    if small_diag is None:
      self._set_max_p(resolution)
      max_p = self.max_p
    else:
      max_p = small_diag[0]


    ax.set_xlim(0,max_p)
    ax.set_ylim(0,max_p)
    self._format_label(ax, resolution, small_diag)
    #ax.invert_yaxis()

  def visualize_single_fig(self, contig, resolution=100000, w_r=8, small_diag=None):
    #figsize = (5,6)
    figsize = (8,9)
    fig = plt.figure(figsize=figsize, facecolor='white', edgecolor='none', frameon=False)
    gs = gridspec.GridSpec(2,1, height_ratios=[w_r,1])
    gs.update(left=0.15, bottom=0.10, top=0.95, right=0.95, wspace=0, hspace=0.01)

    ax = fig.add_subplot(gs[0,0])
    color_ax = fig.add_subplot(gs[-1,:])
    color_ax.set_frame_on(False)

    self.visualize(contig, [ax,], resolution, small_diag=small_diag, color_ax=color_ax, single_color=True, avg_flag=True)
    
    if small_diag is None:
      self._set_max_p(resolution)
      max_p = self.max_p
    else:
      max_p = small_diag[0]
  
    self._format_label(ax, resolution, small_diag)
    ax.set_xlim(0,max_p)
    ax.set_ylim(0,max_p)
    #plt.tight_layout()
    #ax.invert_yaxis()

  def visualize_cuts(self, cuts_fpath, contig, resolution=100000, w_r=5, cut_norm_flag=False, separate_flag=False):
    self._set_max_p(resolution)
    
    with open(cuts_fpath, 'rb') as f:
      cuts = [float(i.strip())/resolution for i in f]
   
    if cut_norm_flag:
      ax_bins = na.arange(self.max_p+1)
      cuts_hist, ax_nbins = na.histogram(cuts,bins=ax_bins)
      cuts_norm = na.outer(cuts_hist, cuts_hist)
    else:
      cuts_norm = None

    bins = 150
    if separate_flag:
      gs = gridspec.GridSpec(4,3, width_ratios=[1,w_r,w_r], height_ratios=[1,w_r,w_r,1])
    else:
      gs = gridspec.GridSpec(3,2, width_ratios=[1,w_r], height_ratios=[1,w_r,1])
    gs.update(left=0.05, bottom=0.05, top=0.95, right=0.95, wspace=0.01, hspace=0.01)

    ax = plt.subplot(gs[1,1])
    histXAx = plt.subplot(gs[0,1], sharex=ax)
    histYAx = plt.subplot(gs[1,0], sharey=ax)
    color_ax = plt.subplot(gs[-1,:])

    if separate_flag:
      ax2 = plt.subplot(gs[1,-1], sharex=ax,sharey=ax)
      ax3 = plt.subplot(gs[-2,1], sharex=ax,sharey=ax)
      self.visualize(contig, [ax,ax2,ax3], resolution, norm_const=cuts_norm, color_ax=color_ax)
    else:
      self.visualize(contig, [ax,], resolution, norm_const=cuts_norm, color_ax=color_ax)
    histXAx.hist(cuts, bins=bins, range=(0,self.max_p), linewidth=0)
    histYAx.hist(cuts, bins=bins, range=(0,self.max_p), orientation='horizontal', linewidth=0)
    
    ax.set_xlim(0,self.max_p)
    ax.set_ylim(0,self.max_p)
    histXAx.set_xticks([])
    histYAx.set_yticks([])
    ax.invert_yaxis()
    histYAx.invert_xaxis()
    
def vis_hic_cut(hic_dat=None):
  contig = "chr19"
  #hic_dat = r"D:\Projects\assembly\data\chr20.dat"
  #cuts_fpath = r"D:\Projects\assembly\data\chr20.hg18.hind3.txt"
  #hic_dat = '/media/T01/data/gm12878/hic/dat/%s.dat'%contig
  if hic_dat is None:
    hic_dat = '/media/T01/data/gm12878/hic/dat/%s_004.dat'%contig
    
  cuts_fpath = "/home/anand/Projects/assembly/data/{0}.hg18.hind3.txt".format(contig)
  
  dat_p = VisualizeDat(hic_dat)
  dat_p.visualize_cuts(cuts_fpath, contig, resolution=100000, cut_norm_flag=False, separate_flag=True)
  plt.show()

def vis_hic(hic_dat, out_fpath=None):
  contig = "chr20"
  dat_p = VisualizeDat(hic_dat)
  dat_p.visualize_single_fig(contig, resolution=100000)
  #dat_p.visualize_linear_fig(contig, resolution=100000)
  if out_fpath is None:
    plt.show()
  else:
    plt.savefig(out_fpath)

def vis_hic_zoom(hic_dat, out_fpath=None, resolution=200, window_size=200, windows=50):
  contig = "chr20"
  #hic_dat = r"D:\Projects\assembly\data\chr20.dat"
  #cuts_fpath = r"D:\Projects\assembly\data\chr20.hg18.hind3.txt"
 
  dat_p = VisualizeDat(hic_dat)
  #dat_p.visualize_single_fig(contig, resolution=resolution, small_diag=(window_size, windows))
  dat_p.visualize_linear_fig(contig, resolution=resolution, small_diag=(window_size, windows))
  
  if out_fpath is None:
    plt.show()
  else:
    plt.savefig(out_fpath)

if __name__ == '__main__':
  import sys
  #vis_hic(True)
  if len(sys.argv)>1 and sys.argv[1].endswith('.dat'):
    try:
      out_fpath = sys.argv[2]
    except IndexError:
      out_fpath = None
    #vis_hic(sys.argv[1], out_fpath)
    #vis_hic_cut(sys.argv[1])
    vis_hic_zoom(sys.argv[1], resolution=200, windows=100, window_size=200)
    #vis_hic_zoom(sys.argv[1], resolution=20, windows=50, window_size=100)
