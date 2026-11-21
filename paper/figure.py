"""Render the manuscript figure from the archived Windows output."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, FixedFormatter

HERE = Path(__file__).resolve().parent
with (HERE.parent/'results'/'diffusion.csv').open() as f:
    rows = [r for r in csv.DictReader(f) if r['density']=='18' and r['vertical_diffusivity']=='1'
            and float(r['width']) < 1 and float(r['lag']) >= 1/512]
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':8.3, 'pdf.fonttype':42})
fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.6, 2.75), width_ratios=[1, 1.85])
fig.subplots_adjust(left=.06, right=.99, top=.86, bottom=.23, wspace=.42)
ax.set_aspect('equal')
ax.add_patch(Rectangle((0, 0), 3, 3, facecolor='#f4f4f1', edgecolor='#757875'))
ax.add_patch(Rectangle((1, 1), 1, 1, facecolor='#cbdce4', edgecolor='#38535e'))
ax.text(1.5, 1.5, 'Hub', ha='center', va='center')
for x,y,label in [(2/3,4/3,'W'),(2,4/3,'E'),(4/3,2/3,'S'),(4/3,2,'N')]:
    ax.add_patch(Rectangle((x, y), 1/3, 1/3, facecolor='#286567'))
    ax.text(x+1/6,y+1/6,label,ha='center',va='center', color='white', fontsize=7)
ax.set(xlim=(-.03,3.03), ylim=(-.03,3.03), xticks=[0,1,2,3], yticks=[0,1,2,3])
ax.set_title('A   Narrow ports', loc='left', fontsize=9, weight='bold')
bx.set_title('B   Cross-port transfer', loc='left', fontsize=9, weight='bold')
t = [float(r['lag']) for r in rows]
for name,label,color,line,marker in [
    ('observed','Observed','#151d24','-','o'), ('bound1','1-state floor','#1a696e','--','s'),
    ('bound3','3-state floor','#a3541d',':','^')]:
    bx.plot(t,[float(r[name]) for r in rows],color=color,ls=line,marker=marker,
            label=label, lw=1.4, markersize=3.5)
bx.set(xscale='log', yscale='log', xlabel='Half-lag  t', ylabel='Stationary transfer mass')
bx.xaxis.set_major_locator(FixedLocator(t))
bx.xaxis.set_major_formatter(FixedFormatter(['1/512','1/128','1/32','1/8','1/2','2']))
bx.minorticks_off()
bx.set_yticks([1e-9,1e-7,1e-5,1e-3,1e-2])
bx.grid(axis='y', color='#dedede', lw=.55)
bx.legend(loc='lower right', frameon=False,fontsize=7.5)
for a in [ax,bx]:
    a.spines[['right','top']].set_visible(False)
    a.tick_params(length=3)
fig.savefig(HERE/'junction.pdf', metadata={'CreationDate':None,'ModDate':None})
