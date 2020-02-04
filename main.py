####
# Visualize Single Cell data
# @author: Zhuoqing Fang
# @email: maxzqfang@stanford.edu
# @version: 0.1
# @time: 2019-12-21
#####
from functools import lru_cache
from os.path import dirname, join
import numpy as np
import pandas as pd
import scanpy as sc 
import colorcet as cc 
from scipy.stats.kde import gaussian_kde
from bokeh.io import curdoc, show
from bokeh.layouts import row, column
from bokeh.models import ColumnDataSource, ColorBar, LinearColorMapper
from bokeh.models import FixedTicker, PrintfTickFormatter
from bokeh.models.glyphs import Patches
from bokeh.models.widgets import Select, TextInput, Dropdown, AutocompleteInput, Div
from bokeh.plotting import figure
from bokeh.palettes import Category10, Spectral6
from bokeh.transform import factor_cmap, factor_mark, linear_cmap, jitter
from bokeh.themes import built_in_themes



DATA_DIR = dirname(__file__)

@lru_cache()
def load_h5ad():
    fname = join(DATA_DIR, "ARPKD.PAGA.intergrated.C123.final.20200108.h5ad")
    #fname = join(DATA_DIR, "old.test.h5ad")
    return sc.read_h5ad(fname)

@lru_cache()
def get_umi(anndat, features):
    #umi = anndat[:, features].X.tolist()
    umi = anndat.obs_vector(features)
    return umi

## Catogorical colors

# https://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
# update 1
# orig reference http://epub.wu.ac.at/1692/1/document.pdf
zeileis_28 = [
    "#023fa5", "#7d87b9", "#bec1d4", "#d6bcc0", "#bb7784", "#8e063b", "#4a6fe3",
    "#8595e1", "#b5bbe3", "#e6afb9", "#e07b91", "#d33f6a", "#11c638", "#8dd593",
    "#c6dec7", "#ead3c6", "#f0b98d", "#ef9708", "#0fcfc0", "#9cded6", "#d5eae7",
    "#f3e1eb", "#f6c4e1", "#f79cd4",
    '#7f7f7f', "#c7c7c7", "#1CE6FF", "#336600",]  # these last ones were added
 
# load data
anndat = load_h5ad() 
#pca = anndat.obsm['X_pca']
pca = anndat.obsm['X_draw_graph_fa']
tsne = anndat.obsm['X_tsne']
umap = anndat.obsm['X_umap']

# set up widgets
#stats = PreText(text='', width=500)
symbol = AutocompleteInput(completions=anndat.var_names.tolist(), 
                           title="Enter Gene Name (e.g. POU5F1 ):", value="AFP")
select = Select(title="Legend:", value="stage_group", options=["clusters","stage","stage_group","donor","Donor","group","res.0.6"])#options=anndat.obs_keys())
# message box
message = Div(text="""Input Gene Name:\n Legend Option: """, width=200, height=100)

## setup data
dd = select.value
umis = get_umi(anndat, symbol.value)
source = ColumnDataSource(data=dict(tSNE1=tsne[:,0].tolist(), 
                                    tSNE2=tsne[:,1].tolist(), 
                                    color=[0]*tsne.shape[0],
                                    umis=[0]*tsne.shape[0],
                                    UMAP1=umap[:,0].tolist(), 
                                    UMAP2=umap[:,1].tolist(),
                                    FA1=pca[:,0].tolist(),
                                    FA2=pca[:,1].tolist(),))
                                    #PC3=pca[:,2].tolist(),))
source_vln = ColumnDataSource(data=dict(xs=[],ys=[], xj=[], yj=[], color=[]))
## setup figures
tools = 'reset,pan,wheel_zoom,box_select,save'
# color_palette= godsnot_102
color_palette = cc.b_glasbey_bw_minc_20
###
catogories = sorted(anndat.obs[dd].unique().tolist())
palette = color_palette[:len(catogories)]
low, high = 0, 1
## transforms
fcmap = factor_cmap('color', palette=palette, factors=[str(c) for c in catogories])
mapper = linear_cmap(field_name='umis', palette="Viridis256",
                     low=low, high=high)
color_bar = ColorBar(color_mapper=mapper['transform'],  
                     width=10, 
                     major_label_text_font_size="10pt",
                     location=(0, 0))

#figures
t1 = figure(plot_width=500, plot_height=500, title="t-SNE", tools=tools)
t1.toolbar.logo = None
t1.xaxis.axis_label = "tSNE1"
t1.yaxis.axis_label = "tSNE2"

#
t2 = figure(plot_width=550, plot_height=500, tools=tools)
t2.toolbar.logo = None
t2.xaxis.axis_label = "tSNE1"
t2.yaxis.axis_label = "tSNE2"

u1 = figure(plot_width=500, plot_height=500, title="UMAP", tools=tools)
u1.toolbar.logo = None
u1.xaxis.axis_label = "UMAP1"
u1.yaxis.axis_label = "UMAP2"

u2 = figure(plot_width=550, plot_height=500, tools=tools)
u2.toolbar.logo = None
u2.xaxis.axis_label = "UMAP1"
u2.yaxis.axis_label = "UMAP2"

## trajectory
p1 = figure(plot_width=500, plot_height=500, title="Trajectory Inference", tools=tools)
p1.toolbar.logo = None
p1.xaxis.axis_label = "FA1"
p1.yaxis.axis_label = "FA2"

p2 = figure(plot_width=550, plot_height=500, tools=tools)
p2.toolbar.logo = None
p2.xaxis.axis_label = "FA1"
p2.yaxis.axis_label = "FA2"

########### clustering plots ##################
## tSNE
t1.circle('tSNE1', 'tSNE2', size=3, source=source, legend_field="color", color= fcmap,
            selection_color="orange", alpha=0.6, nonselection_alpha=0.1, selection_alpha=0.4)
## UMAP
u1.circle('UMAP1', 'UMAP2', size=3, source=source, legend_field="color", color= fcmap,
            selection_color="orange", alpha=0.6, nonselection_alpha=0.1, selection_alpha=0.4)
## PCA
p1.circle('FA1', 'FA2', size=3, source=source, legend_field="color", color= fcmap,
            selection_color="orange", alpha=0.6, nonselection_alpha=0.1, selection_alpha=0.4)   

###### gene expression plots ##################
# tSNE
t2.scatter('tSNE1', 'tSNE2', size=3, source=source, 
            fill_color=mapper,
            line_color=None, selection_color="orange", alpha=0.6, 
            nonselection_alpha=0.1, selection_alpha=0.4)
t2.add_layout(color_bar, 'right')

# UMAP
u2.scatter('UMAP1', 'UMAP2', size=3, source=source, 
            fill_color=mapper,
            line_color=None, selection_color="orange", alpha=0.6, 
            nonselection_alpha=0.1, selection_alpha=0.4)
u2.add_layout(color_bar, 'right') 
#PCA
p2.scatter('FA1', 'FA2', size=3, source=source, 
            fill_color=mapper,
            line_color=None, selection_color="orange", 
            alpha=0.6, nonselection_alpha=0.1, selection_alpha=0.4)
p2.add_layout(color_bar, 'right')


## voilinplot
volin = figure(plot_width=1000, plot_height=500, 
               tools = 'reset, pan,wheel_zoom, save')
volin.patches(xs='xs', ys='ys', alpha=0.6, fill_color='color', line_color='black', source=source_vln)
#volin.circle(x=jitter('xj', 0.4), y='yj', size=2, color='black', alpha=0.4, source=source_vln)

volin.toolbar.logo = None
volin.yaxis.axis_label = "Expression"

volin.outline_line_color = None
volin.background_fill_color = "#efefef"
volin.xaxis.major_label_orientation = np.pi/4
volin.xgrid.grid_line_color = None
volin.ygrid.grid_line_color = "#dddddd"
volin.axis.minor_tick_line_color = None
volin.axis.major_tick_line_color = None
volin.axis.axis_line_color = None

def volin_change(gene, catogory, umis, bins=1000, cut=2):
    # update axis
    cats = sorted(catogory.unique())
    x_range = (np.arange(0, len(cats)) * 3).tolist()
    ## update data
    color = color_palette[:len(cats)]
    xs = []
    ys =  []
    xj = []
    yj = []

    for i, cat in zip(x_range, cats):
        umi = umis[catogory == cat]
        kde = gaussian_kde(umi) 
        # same default paramter same as seaborn.violinplot
        bw_used = kde.factor * umi.std(ddof=1) * cut
        support_min = umi.min() - bw_used 
        support_max = umi.max() + bw_used
        kde_support = np.linspace(support_min, support_max, bins)

        # generate y and x values
        yy = np.concatenate([kde_support, kde_support[::-1]])  # note: think about how bokeh patches work 
        x = kde(kde_support)  ### you may change this if the violin not fit in
        x /= x.max() # scale the relative area under the kde curve, resulting max density will be 1
        x2 = -x
        xx = np.concatenate([x, x2[::-1]]) + i
        xs.append(xx)
        ys.append(yy)
        #xj.append([i]*len(umi))
        #yj.append(umi)
      
    source_vln.data = dict(xs=xs, ys=ys, color=color)
    #source_vln.data = dict(xs=xs, ys=ys, color=color, xj=xj, yj=yj)

    volin.xaxis.ticker = FixedTicker(ticks= x_range)
    volin.xaxis.major_label_overrides = {k: str(v) for k, v in zip(x_range, cats)}
    volin.title.text = gene

# set up callbacks
def gene_change():
    ## update text input and select attribute
    gene = symbol.value.strip()
    dd = select.value

    message.text = "Input Gene Name: {g} \nLegend Option: {cat}".format(g=gene,cat=dd)
    # select gene expression value
    umis = get_umi(anndat, gene)
    ## update existing tranform
    vmin, vmax = np.percentile(umis, [2, 98])
    mapper['transform'].low= vmin
    mapper['transform'].high= vmax
    ## update title
    p2.title.text = gene 
    u2.title.text = gene
    t2.title.text = gene   
    ## update source data
    source.data.update(umis=umis,)
    
    # update volin
    clusters = anndat.obs[dd]
    volin_change(gene, clusters, umis, bins=1000)


# set up callbacks
def factor_change():
    ## update text input and select attribute
    gene = symbol.value.strip()
    dd = select.value
 
    message.text = "Input Gene Name: {g} \nLegend Option: {cat}".format(g=gene,cat=dd)
    # select gene expression value
    umis = get_umi(anndat, gene) 
    # update factor color 
    clusters = anndat.obs[dd]
    cats = sorted(clusters.unique().tolist())
    fcmap['transform'].factors= [ str(c) for c in cats] 
    fcmap['transform'].palette=color_palette[:len(cats)]
    ## update source data
    source.data.update(color=clusters.astype(str).tolist())
    ## update violin
    volin_change(gene, clusters, umis, bins=1000)

### input on change
symbol.on_change('value', lambda attr, old, new: gene_change())
select.on_change('value', lambda attr, old, new: factor_change())

# set up layout
sb = column(symbol, select, message)

series_t = row(t1, t2, sb)
series_u = row(u1, u2, volin)
series_p = row(p1, p2, )
layout = column(series_t, series_u, series_p,)

# initialize
#update()
gene_change()
factor_change()

curdoc().add_root(layout)
#curdoc.theme = 'dark_minimal'
curdoc().title = "ARPKD"
