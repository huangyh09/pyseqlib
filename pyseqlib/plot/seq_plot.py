# This file is to produce figures for gene struncture and RNA-seq reads
# density. It is largely borrowed from Yarden's sashimi_plot

import numpy as np
import matplotlib.pyplot as plt

def set_frame(ax):
    """Setting the frame of the figure box for plot gene.
    """
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    return ax

def plot_gene(g, tran_colors=[], exon_emph=[], exon_colors=[], 
              gap=2.0, narrows=40, nxticks=5):
    """
    Draw the gene structure. 
    Note: exonw_width is 1.0, y_loc_base is 0.
    
    Parameters
    ----------
    g: 
        gene in brie package
    tran_colors: 
        a list of colors for transcripts
    exon_emph:
        a list of exons, e.g., [[4,8], [12,20]]
    exon_colors:
        a list of colors for exons to emphasize
    gap: 
        float, the gap between two transcripts (width is 1)
    narrows:
        int, the number of arrows
    nxticks:
        int, the number of xticks
    """
    mRNAs = [t.exons for t in g.trans]
    for i in range(len(mRNAs)):
        t_color = "k"
        if len(tran_colors) > i: 
            t_color = tran_colors[i]
        mRNA = mRNAs[i]
        yloc = gap*i
        
        # plot exon
        for j in range(len(mRNA)):
            x = [mRNA[j][0], mRNA[j][1], mRNA[j][1], mRNA[j][0]]
            y = [yloc-0.5, yloc-0.5, yloc+0.5, yloc+0.5]
            
            _color = t_color
            if exon_emph.count([mRNA[j][0], mRNA[j][1]]) == 1:
                idx = exon_emph.index([mRNA[j][0], mRNA[j][1]])
                _color = exon_colors[idx]
            plt.fill(x, y, color=_color, edgecolor=_color, zorder=20)
                
        # plot intron.
        xmin = np.min(mRNA)
        xmax = np.max(mRNA)
        plt.plot([xmin, xmax], [yloc, yloc], color='k', lw=.5)

        # plot intron arrows.
        spread = 0.2 * (g.stop-g.start) / (narrows+1)
        xloc = np.linspace(g.start, g.stop, narrows)
        for i in range(narrows):
            if xloc[i] < xmin or xloc[i] > xmax:
                continue  
            if g.strand == '+':
                x = [xloc[i]-spread, xloc[i], xloc[i]-spread]
            else:
                x = [xloc[i]+spread, xloc[i], xloc[i]+spread]
            y = [yloc-0.2, yloc, yloc+0.2]
            plt.plot(x, y, lw=.5, color='k') 
            
    #plot chromosome axis
    plt.yticks([])
    xticks = np.linspace(g.start, g.stop, nxticks).astype(int)
    plt.xticks(xticks, xticks)
    plt.xlabel('Genomic coordinate (%s), "%s" strand'%(g.chrom, g.strand))
    
    xlim1 = g.start - (g.stop - g.start) * 0.025
    xlim2 = g.stop  + (g.stop - g.start) * 0.025
    plt.xlim(xlim1, xlim2) 
    plt.ylim(-gap+0.5, len(mRNAs)*gap-0.5)









# def plot_gene(tx_start, mRNAs, strand,  graphcoords, reverse_minus=False, 
#               exonwidth=0.3, narrows=50, gap=1.5, bound=0.0, ybase=0):
#     """
#     Draw the gene structure.
#     """
#     cnt = 0
#     for i in range(len(mRNAs)):
#         mRNA = mRNAs[i]
#         yloc = gap*exonwidth*cnt + ybase
#         if len(mRNA) == 2: continue
#         cnt += 1
#         exon_sorted = np.sort(mRNA, axis=0)
#         for j in range(len(mRNA)):
#             s = mRNA[j][0] - tx_start
#             e = mRNA[j][1] - tx_start
#             x = [graphcoords[s], graphcoords[e], graphcoords[e], graphcoords[s]]
#             y = [yloc - exonwidth/2, yloc - exonwidth/2,
#                  yloc + exonwidth/2, yloc + exonwidth/2]
#             fill(x, y, 'k', lw=.5, zorder=20)
#             if len(mRNA) == 3 and mRNA[j][0] == exon_sorted[1,0]: 
#                 fill(x, y, 'gray', lw=.0, zorder=20)
#             else:
#                 fill(x, y, 'k', lw=.5, zorder=20)

#         # Draw intron.
#         axhline(yloc, color='k', lw=.5)

#         # Draw intron arrows.
#         spread = .2 * max(graphcoords) / narrows
#         for i in range(narrows):
#             loc = float(i) * max(graphcoords) / narrows
#             if strand == '+' or reverse_minus:
#                 x = [loc-spread, loc, loc-spread]
#             else:
#                 x = [loc+spread, loc, loc+spread]
#             y = [yloc - exonwidth/5, yloc, yloc + exonwidth/5]
#             plot(x, y, lw=.5, color='k')