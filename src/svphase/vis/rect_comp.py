import numpy as na
from itertools import cycle

from matplotlib import pyplot as plt

import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection


colors = {'b': '#0000CD', 'both': '#9400D3', 'a': '#B22222', 'mis':0.4, 'cor':'0.7'}

class Rect(object):
  def __init__(self, px, py, w, h):
    self.px = px
    self.py = py
    self.w = w
    self.h = h
 
class RectVertiError(Rect):
  def __init__(self, px, py, w, h):
    Rect.__init__(self, px, py, w, h)
    self.labels = ['Correct',]
    self.colors = map(colors.__getitem__, ['cor', ])  
  def get_patches(self, vconts):
    return [mpatches.Rectangle((self.px,py),self.w,h, ec='none', fc=self.colors[idx], alpha=0.9) for idx, (py, h) in enumerate(vconts)]
  def draw(self, class_s, tot):
    #tot = float(class_s + error_s)
    class_h = self.h*class_s/float(tot)
    #error_h = self.h*error_s/tot
    return [(self.py, class_h),]

class RectVertiComp(RectVertiError):
  def __init__(self, px, py, w, h, labels):
    RectVertiError.__init__(self, px, py, w, h)
    #self.labels = ['HiC-only', 'Both', 'WGS-only']
    self.labels = labels
    self.colors = map(colors.__getitem__, ['a', 'both', 'b'])
  def draw(self, a_s, both_s, b_s):
    tot = float(a_s + b_s + both_s)
    a_h = self.h*a_s/tot
    both_h = self.h*both_s/tot
    b_h = self.h*b_s/tot

    return [(self.py, a_h),(self.py + a_h, both_h), (self.py+a_h+both_h, b_h)]

class RectHorizComp(Rect):
  def __init__(self,px,py,w,h):
    Rect.__init__(self, px, py, w, h)
    self.cont = []
    self.labels=["pA", "pB", "homA", "nan", "inc"]
    #self.colors = ['%.2f'%i for i in na.linspace(0.10,0.5,len(self.labels))]
    self.colors = ['0.9']*len(self.labels)
  def get_patches(self, hconts):
    return [mpatches.Rectangle((px,self.py),w, self.h, ec='none', fc=self.colors[idx])  for idx, (px, w) in enumerate(hconts) if w>0]
  def draw(self, phased_a_s, phased_b_s, homo_s, nan_s, incor_s, wspace_frac=0):
    tot = float(phased_a_s + phased_b_s + homo_s + nan_s + incor_s)    
    wspace_frac *= self.w
    phased_a_w = self.w*phased_a_s/tot
    phased_b_w = self.w*phased_b_s/tot
    homo_w = self.w*homo_s/tot
    nan_w = self.w*nan_s/tot
    incor_w = self.w*incor_s/tot

    x = [self.px, phased_a_w, phased_b_w, homo_w, nan_w,incor_w]
    # returns containers
    return zip(na.cumsum(x[:-1]),na.array(x[1:])-wspace_frac)

class RectHorizPhaseMergedComp(RectHorizComp):
  def __init__(self,px,py,w,h):
    RectHorizComp.__init__(self, px, py, w, h)
    self.cont = []
    self.labels=["Phased", "Homozygous", "No data", "Inc."]
    self.colors = ['%.2f'%i for i in na.linspace(0.10,0.5,4)]
  def draw(self, phased_s, homo_s, nan_s, incor_s, wspace_frac=0):
    tot = float(phased_s + homo_s + nan_s + incor_s)    
    wspace_frac *= self.w
    phased_w = self.w*phased_s/tot
    homo_w = self.w*homo_s/tot
    nan_w = self.w*nan_s/tot
    incor_w = self.w*incor_s/tot
    # returns containers
    return [(self.px, phased_w-wspace_frac),
            (self.px+phased_w, homo_w-wspace_frac),
            (self.px+phased_w+homo_w, nan_w-wspace_frac), 
            (self.px+phased_w+homo_w+nan_w, incor_w-wspace_frac)]

def model_diff_ax(ax, c,s, text_ylabels):
  patches = []
  horiz_cont = RectHorizComp(0,0, 2,1)

  try:
    nan_s = s.ix['nan']
  except KeyError:
    nan_s = 0
  try:
    homo_s = s.ix['homA']
  except KeyError:
    homo_s = 0
  try:
    incor_s = s.ix['inc']
  except KeyError:
    incor_s = 0

  hconts = horiz_cont.draw(phased_a_s=s.ix['pA'],phased_b_s=s.ix['pB'],
                           homo_s=homo_s, nan_s=0, incor_s=incor_s,wspace_frac=0.01) 
  patches = horiz_cont.get_patches(hconts)

  vertis = []
  for lbl, (px, w) in zip(horiz_cont.labels, hconts):
    if w<0: continue
    v = RectVertiComp(px, horiz_cont.py, w, horiz_cont.h, labels=text_ylabels)
    vconts = v.draw(a_s=c['a'].ix[lbl], both_s=c['both'].ix[lbl], b_s=c['b'].ix[lbl])
    patches.extend(v.get_patches(vconts))
    vertis.append((v,vconts))
  #print "Vertical Containers"
  #print vertis
  color = '0.1'

  #collection = PatchCollection(patches)
  for p in patches:
    ax.add_patch(p)
  #ax.add_collection(collection)
  lbl_fix = dict(zip(s.index, s.index))
  lbl_fix['homA'] = 'hom.'
  lbl_fix['inc'] = 'incorrect'
  for (px,w), lbl in zip(hconts,horiz_cont.labels):
    if w<0: continue
    #print px, w, lbl, color
    ax.text(px + w/2, 0, lbl_fix[lbl], color=color, ha='center', va='baseline', size='large')
  for (py, h), lbl, c in zip(vertis[0][1],vertis[0][0].labels, vertis[0][0].colors):
    ax.text(2, py+h/2 , lbl, ha='left', va='center', color=c, size='large')
  
  ax.set_ylabel('Fraction of correct predictions')
  ax.set_xlim(horiz_cont.px,horiz_cont.px+horiz_cont.w)
  ax.set_ylim(horiz_cont.py,horiz_cont.py+horiz_cont.h)

  ax.set_frame_on(False)
  ax.get_yaxis().tick_left()
  #ax.axes.get_yaxis().set_visible(False)
  ax.axes.get_xaxis().set_visible(False)

def error_ax(ax, c, s, tick_flag=True, min_accuracy=0.1, background_flag=True): 
  patches = []
  horiz_cont = RectHorizComp(0,0, 2,1)
  try:
    nan_s = s.ix['nan']
  except KeyError:
    nan_s = 0
  try:
    homo_s = s.ix['homA']
  except KeyError:
    homo_s = 0
  try:
    incor_s = s.ix['inc']
  except KeyError:
    incor_s = 0

  hconts = horiz_cont.draw(phased_a_s=s.ix['pA'],phased_b_s=s.ix['pB'],
                           homo_s=homo_s, nan_s=nan_s, incor_s=incor_s,wspace_frac=0.01) 
  if background_flag:
    patches = horiz_cont.get_patches(hconts)
  else:
    patches = []

  vertis = []
  for lbl, (px, w) in zip(horiz_cont.labels, hconts):
    if w<0: continue
    v = RectVertiError(px, horiz_cont.py, w, horiz_cont.h)
    vconts = v.draw(class_s=c[lbl], tot=s[lbl])
    patches.extend(v.get_patches(vconts))
    vertis.append((v,vconts))
  #print "Vertical Containers"
  #print vertis

  #collection = PatchCollection(patches)
  for p in patches:
    ax.add_patch(p)
  #ax.add_collection(collection)
  #for (px,w), lbl, color in zip(hconts,horiz_cont.labels, horiz_cont.colors):
  color = '0.1'
  
  lbl_fix = dict(zip(s.index, s.index))
  lbl_fix['homA'] = 'hom.'
  lbl_fix['inc'] = 'incorrect'
  lbl_fix['nan'] = 'excl.'

  for (px,w), lbl in zip(hconts,horiz_cont.labels):
    if w<0: continue
    #print px, w, lbl, color
    #ax.text(px + w/2, 1, lbl, color=color, ha='center')
    cfrac = c[lbl]/float(s[lbl])-0.05
    if lbl!='nan':
      if cfrac>min_accuracy:
        ax.text(px + w/2, cfrac, "{0:d}/{1:d}".format(c[lbl],s[lbl]), color=color, ha='center', va='top', size='medium')
      ax.text(px + w/2, min_accuracy, lbl_fix[lbl], color=color, ha='center', va='bottom', size='large')
    else:
      ax.text(px + w/2, min_accuracy, "{0}\n{1}".format(s[lbl],lbl_fix[lbl]), color=color, ha='center', va='bottom')
  #for (py, h), lbl, c in zip(vertis[0][1],vertis[0][0].labels, vertis[0][0].colors):
  #  ax.text(0, py+h/2 , lbl, ha='right', va='center', color=c)

  ax.set_xlim(horiz_cont.px,horiz_cont.px+horiz_cont.w)
  ax.set_ylim(horiz_cont.py+min_accuracy,horiz_cont.py+horiz_cont.h)

  ax.set_frame_on(False)
  if tick_flag:
    ax.get_yaxis().tick_left()
  ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%0.2f'))
  ax.axes.get_xaxis().set_visible(False)
 
def test():
  patches = []
  horiz_cont = RectPhaseMergedHorizComp(0,0, 2,1)
  hconts = horiz_cont.draw(30,10,5,10,wspace_frac=0.01) 
  patches.extend([mpatches.Rectangle((px,horiz_cont.py),w, horiz_cont.h, ec='none', fc=horiz_cont.colors[idx])  for idx, (px, w) in enumerate(hconts)])
  #print "Horizontal Containers"
  #print hconts

  ratios = cycle([(40,30,10), (30,50,5)])
  vertis = []
  for lbl, (px, w) in zip(horiz_cont.labels, hconts):
    if lbl=='No data': continue
    v = RectVertiComp(px, horiz_cont.py, w, horiz_cont.h)
    hic,both,wgs = ratios.next()
    vconts = v.draw(hic, both, wgs)
    patches.extend([mpatches.Rectangle((v.px,py),v.w,h, ec='none', fc=v.colors[idx], alpha=0.5) for idx, (py, h) in enumerate(vconts)])
    vertis.append((v,vconts))
  #print "Vertical Containers"
  #print vertis


  fig = plt.figure(facecolor='white', edgecolor='none')

  ax = fig.add_subplot(111)
  #collection = PatchCollection(patches)
  for p in patches:
    ax.add_patch(p)
  #ax.add_collection(collection)
  for (px,w), lbl, c in zip(hconts,horiz_cont.labels, horiz_cont.colors):
    #print px, w, lbl, c
    ax.text(px + w/2, 1, lbl, color=c, ha='center')
  for (py, h), lbl, c in zip(vertis[0][1],vertis[0][0].labels, vertis[0][0].colors):
    ax.text(0, py+h/2 , lbl, ha='right', va='center', color=c)

  plt.xlim(horiz_cont.px,horiz_cont.px+horiz_cont.w)
  plt.ylim(horiz_cont.py,horiz_cont.py+horiz_cont.h)
  #plt.axis('equal')
  plt.axis('off')
  plt.show()
if __name__ == '__main__':
  test() 
