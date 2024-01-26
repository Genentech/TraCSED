#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#Temporal_Features_30x


# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler

import os
#var = 'barmodels_pcares0.25_BREAST75'
#os.chdir('/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/'+var)
#os.listdir()


# In[2]:


def tempFeaturesPerBar(barcode):
    df = pd.read_csv(barcode)
    basedf = df[df['feature'] == 'BASELINE']
    realp = eval(basedf['prediction'].tolist()[0])
    realt = eval(basedf['true'].tolist()[0])
    featerr = {}
    for feature in df['feature'].tolist():
        if feature != 'BASELINE':
            smalldf = df[df['feature'] == feature]
            pdf = pd.DataFrame()
            for col in smalldf.columns.tolist():
                if 'prediction_' in col:
                    pdf[col] = eval(smalldf[col].tolist()[0])
        
            pdf['Real_Prediction'] = realp
            pdf['Real_True'] = realt
            pdf2 = pdf
            pdf2 = pdf2.sub(pdf2['Real_Prediction'], axis=0)

            predictions = []
            #true = realt
            #realpred = realp
            for each in pdf2.index.tolist():
                predictions.append(np.mean(pdf2.iloc[each].tolist()[:30]))
            err = []
            for i in range(len(predictions)):
                if realt[i] > predictions[i] and predictions[i] > realp[i]:
                    err.append(-abs(predictions[i])) # permutation is closer to true than model, prediction is off low
                if realt[i] < predictions[i] and predictions[i] < realp[i]:
                    err.append(-abs(predictions[i])) # permutation is closer to true than model, prediction is off high
                if realt[i] > predictions[i] and predictions[i] < realp[i]:
                    err.append(abs(predictions[i]))
                if realt[i] < predictions[i] and predictions[i] > realp[i]:
                    err.append(abs(predictions[i]))
            featerr[feature] = err
            
    return featerr


# In[3]:


var = 'barmodels_palbo_plsrres0.3_HALL3'
os.chdir('/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/'+var)
barcode = 'transformer_fullGFPBC_libB_92196_feature_importance_full.csv'
feature = 'Score.HALLMARK_E2F_TARGETS.1'


df = pd.read_csv(barcode)
basedf = df[df['feature'] == 'BASELINE']
realp = eval(basedf['prediction'].tolist()[0])
realt = eval(basedf['true'].tolist()[0])
featerr = {}

smalldf = df[df['feature'] == feature]
pdf = pd.DataFrame()
for col in smalldf.columns.tolist():
    if 'prediction_' in col:
        pdf[col] = eval(smalldf[col].tolist()[0])

pdf['Real_Prediction'] = realp
pdf['Real_True'] = realt

predictions = []
for each in pdf.index.tolist():
    predictions.append(np.mean(pdf.iloc[each].tolist()[:30]))
    
plt.figure(figsize=(15,5))
plt.plot(range(len(realt)), realt, color='black')
plt.plot(range(len(predictions)), predictions, color='red', alpha=0.5)
plt.show()


# In[4]:


def fullPlot(bardict):
    import seaborn as sns
    get_ipython().run_line_magic('matplotlib', 'inline')

    bardf = pd.DataFrame(bardict)
    c = sns.clustermap(bardf.T, annot=False, xticklabels=True, row_cluster=True, col_cluster=False,
                   yticklabels=True, center=0, cmap='coolwarm')
    c.tick_params(labelsize=6)
    plt.show()
    
    
def prunedPlot(bardict, prune_val):
    import seaborn as sns
    get_ipython().run_line_magic('matplotlib', 'inline')

    bardf = pd.DataFrame(bardict)
    print('before prune', bardf.shape)
    bardf2 = pd.DataFrame()
    for col in bardf.columns.tolist():
        if max(bardf[col]) >= prune_val:
            bardf2[col] = bardf[col]
    print('after prune', bardf2.shape)

    d = sns.clustermap(bardf2.T, annot=False, xticklabels=True, row_cluster=True, col_cluster=False,
                       yticklabels=True, center=0, cmap='coolwarm')
    d.tick_params(labelsize=6)
    plt.show()
    
    
def prunedPlot2(bardict, cutoff):
    import seaborn as sns
    get_ipython().run_line_magic('matplotlib', 'inline')

    bardf = pd.DataFrame(bardict)
    print('before prune', bardf.shape)
    bardf2 = pd.DataFrame()
    mcols = []
    for col in bardf.columns.tolist():
        mcols.append(max(bardf[col]))
    mcols.sort(reverse=True)
    topcut = []
    for col in bardf.columns.tolist():
        if max(bardf[col]) >= mcols[cutoff-1]:
            bardf2[col] = bardf[col]
            topcut.append(col)
    print('after prune', bardf2.shape)

    d = sns.clustermap(bardf2.T, annot=False, xticklabels=True, row_cluster=True, col_cluster=False,
                       yticklabels=True, center=0, cmap='coolwarm')
    d.tick_params(labelsize=6)
    plt.show()
    return topcut


# In[5]:


var = 'barmodels_combo_plsrres0.5_HALL3'
os.chdir('/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/'+var)

alltoppaths = []
for f in os.listdir():
    if '_fullGFPBC_libB_' in f:
        bardict = tempFeaturesPerBar(f)
        print(f)
        #prunedPlot(bardict, 0.04)
        top = prunedPlot2(bardict, 30)
        alltoppaths.append(top)


# In[6]:


countpaths = {}
for l in alltoppaths:
    for p in l:
        if p in countpaths.keys():
            countpaths[p] += 1
        else:
            countpaths[p] = 1
countpathsl = sorted(countpaths.items(), key=lambda x:x[1], reverse=True)
#print(countpathsl)


# In[4]:


var = 'barmodels_9545_plsrres0.3_GENE3'
os.chdir('/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/'+var)

alltopgenes = []
for f in os.listdir():
    if '_fullGFPBC_libB_' in f:
        print(f)
        bardict = tempFeaturesPerBar(f)
        #prunedPlot(bardict, 0.005)
        top = prunedPlot2(bardict, 40)
        alltopgenes.append(top)


# In[5]:


countgenes = {}
for l in alltopgenes:
    for p in l:
        if p in countgenes.keys():
            countgenes[p] += 1
        else:
            countgenes[p] = 1
countgenesl = sorted(countgenes.items(), key=lambda x:x[1], reverse=True)
print(countgenesl[:50])


# In[7]:


os.chdir('/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/')
#jfile = 'tran_9545_pcares0.25_BREAST.json'
jfile = 'json_files/tran_combo_plsrres0.5_HALL3.json'
#os.listdir()

import json
f = open(jfile)
data = json.load(f)

barcodes = [i for i in data.keys()]
new_barcodes = []
for i in barcodes:
    df = pd.DataFrame(data[i])
    #print(i, df.shape[0])
    if df.shape[0] > 200: # 50 or 60
        new_barcodes.append(i)
print(new_barcodes)
print(len(new_barcodes))

bar_dict = {new_barcodes[i]:i for i in range(len(new_barcodes))}


# In[8]:


df = pd.DataFrame(data['GFPBC_libB_90850'])
df.iloc[:, 0:len(df.columns)-1]


# In[9]:


import matplotlib.pyplot as plt

def plotFeaturePerBar(feat, jdata, bars):
    for b in bars:
        df = pd.DataFrame(jdata[b])
        featbar = df[feat].tolist() ######## not normalized
        featbar = [i / max(featbar) for i in featbar]
        ptime = df['Pseudotime'].tolist()
        plt.plot(ptime, featbar, label=b)
        plt.xlabel('Pseudotime')
        plt.ylabel('Max-Normalized Score')
        plt.title('Barcodes Activity of '+feat)
        plt.legend()
    plt.show()


# In[10]:


def plotAvgFeaturePerBar(feat, jdata, bars, seqlen, resbars, timezoom):
    
    c=['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'yellow', 'black', 'yellow']
    colors = {bars[i]:c[i] for i in range(len(bars))}
    
    for b in bars:
        df = pd.DataFrame(jdata[b])
        featbar = df[feat].tolist() ######## not normalized
        # normalize feature or not
        #featbar = [i / max(featbar) for i in featbar]
        ptime = df['Pseudotime'].tolist()
        ##### mean of 8-20 sequence
        seqnormscore = []
        for i in range(seqlen, len(featbar)+seqlen):
            seqnormscore.append(np.mean(featbar[i-seqlen:i]))
        if b in resbars:
            s = 'dotted'
            a = 1
        else:
            s = 'solid'
            a = 0.25
        barname = b.split('_')[-1]
        plt.plot(ptime, seqnormscore, linestyle=s, label=barname, color=colors[b], alpha=a)
        plt.xlabel('Pseudotime')
        plt.ylabel('Max-Normalized Avg of Seq'+str(seqlen))
        plt.title('Barcodes Activity of '+feat)
        plt.xlim(timezoom)
        plt.legend(loc='best', fontsize=6)
    plt.show()


# In[10]:


# top 20 paths
resbars9545 = ['GFPBC_libB_90850', 'GFPBC_libB_10678', 'GFPBC_libB_92196', 'GFPBC_libB_10737']
resbarspalbo = ['GFPBC_libB_90850', 'GFPBC_libB_92196']
resbarscombo = ['GFPBC_libB_92196', 'GFPBC_libB_37798']

for p in countpathsl[:20]:
    plotAvgFeaturePerBar(p[0], data, new_barcodes, 50, resbarscombo, [0,1])


# In[14]:


# survival
#plotAvgFeaturePerBar(var, jdata, barcodes, seqlen_act, resbars, timezoom)
resbars9545 = ['GFPBC_libB_90850', 'GFPBC_libB_10678', 'GFPBC_libB_92196', 'GFPBC_libB_10737']
resbarspalbo = ['GFPBC_libB_90850', 'GFPBC_libB_92196']
resbarscombo = ['GFPBC_libB_92196', 'GFPBC_libB_37798']

plotAvgFeaturePerBar('Survival', data, new_barcodes, 1, resbarscombo, [0,1])


# In[11]:


def plotFeatNormModel(bars, var, seqlen_mod, resbars, master, bar_ptime):
    path_dict = {}
    for b in bars:
        df = pd.DataFrame(master[b])
        df2 = df.T
        df2[df2 < 0] = 0
        scaler = MinMaxScaler(feature_range=(0, 1)) 
        df2 = pd.DataFrame(scaler.fit_transform(df2))
        df2.index = df.columns

        df = df2.T
        path_dict[b] = df[var].tolist()
    
    smooth_dict = {}
    for bar in path_dict.keys():
        featbar = path_dict[bar]
        seqnormscore = []
        for i in range(seqlen_mod, len(featbar)+seqlen_mod):
            seqnormscore.append(np.mean(featbar[i-seqlen_mod:i]))
        smooth_dict[bar] = seqnormscore
        
    for bar in smooth_dict.keys():
        if bar in resbars:
            s = 'dotted'
        else:
            s = 'solid'
        #featbar = [i / max(smooth_dict[bar]) for i in smooth_dict[bar]]
        plt.plot(bar_ptime[bar][8:], smooth_dict[bar], linestyle=s, label=bar)
    plt.xlabel('Pseudotime')
    plt.ylabel('FeatNorm Avg of Seq'+str(seqlen_mod))
    plt.title('FeatNorm Model Importance of '+var)
    plt.legend(loc='best', fontsize=6)
    plt.ylim(bottom=0)
    plt.show()


# In[12]:


import pandas as pd
import numpy as np
from sklearn.neighbors import KernelDensity

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as grid_spec


def ridgePlot(path):

    times = [x for x in np.unique(df.RealTime)]
    days = ['Day 0', 'Day 1', 'Day 4', 'Day 8', 'Day 26']
    #print(times)
    colors = ['#0000ff', '#3300cc', '#660099', '#990066', '#cc0033', '#ff0000']

    gs = grid_spec.GridSpec(len(times),1)
    fig = plt.figure(figsize=(7.5,2.5))

    i = 0

    ax_objs = []
    for t in times:
        t = times[i]
        x = np.array(df[df.RealTime == t].Pseudotime)
        x_d = np.linspace(0,1, 1000) ######

        kde = KernelDensity(bandwidth=0.03, kernel='gaussian')
        kde.fit(x[:, None])

        logprob = kde.score_samples(x_d[:, None])

        # creating new axes object
        ax_objs.append(fig.add_subplot(gs[i:i+1, 0:]))

        # plotting the distribution
        ax_objs[-1].plot(x_d, np.exp(logprob),color="#f0f0f0",lw=1)
        ax_objs[-1].fill_between(x_d, np.exp(logprob), alpha=1,color=colors[i])


        # setting uniform x and y lims
        ax_objs[-1].set_xlim(0,1)
        #ax_objs[-1].set_ylim(0,2.5)

        # make background transparent
        rect = ax_objs[-1].patch
        rect.set_alpha(0)

        # remove borders, axis ticks, and labels
        ax_objs[-1].set_yticklabels([])

        if i == len(times)-1:
            ax_objs[-1].set_xlabel("Pseudotime", fontsize=16,fontweight="bold")
        else:
            ax_objs[-1].set_xticklabels([])

        spines = ["top","right","left","bottom"]
        for s in spines:
            ax_objs[-1].spines[s].set_visible(False)

        #adj_t = t.replace(" ","\n")
        ax_objs[-1].text(-0.02,0,days[i],fontweight="bold",fontsize=14,ha="right")


        i += 1

    gs.update(hspace=-0.7)

    fig.text(0.07,0.85,"Distribution of Experiment along Pseudotime",fontsize=14)

    plt.tight_layout()
    if path != '':
        plt.savefig(path+'RidgePlot.png')
    plt.show()
    return ax_objs


# In[13]:


# pathway model scores across ptime on xaxis for each barcode in heatmap form
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

def masterBardict(barcodes, path):
    
    # get bardict for each barcode
    master_bardict = {}
    for bar in barcodes:
        barstring = path+'transformer_full'+bar+'_feature_importance_full.csv'
        bardict = tempFeaturesPerBar(barstring)
        master_bardict[bar] = bardict
    return master_bardict

def plotModelScoresPerBar(var, barcodes, jdata, seqlen_mod, seqlen_act, master, resbars, timezoom):
    
    c=['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'yellow', 'black', 'yellow']
    colors = {barcodes[i]:c[i] for i in range(len(barcodes))}
    
    # get bardict for each barcode
    path_bardict = {}
    bar_ptime = {}
    for bar in barcodes:
        # lookup and store the error scores for path/gene
        path_bardict[bar] = master[bar][var]
        
        # get the ptime for cells in barcode from jdata
        df = pd.DataFrame(jdata[bar])
        bar_ptime[bar] = df['Pseudotime'].tolist()
    
    # smooth scores along ptime after normalized
    ##### mean of seqlen
    smooth_dict = {}
    for bar in path_bardict.keys():
        featbar = path_bardict[bar]
        seqnormscore = []
        for i in range(seqlen_mod, len(featbar)+seqlen_mod):
            seqnormscore.append(np.mean(featbar[i-seqlen_mod:i]))
        smooth_dict[bar] = seqnormscore
    
    # plot model data as max-norm plot
    #resbars = ['GFPBC_libB_90850', 'GFPBC_libB_10678', 'GFPBC_libB_92196']
    for bar in smooth_dict.keys():
        if bar in resbars:
            s = 'dotted'
            a = 1
        else:
            s = 'solid'
            a = 0.25
        #featbar = [i / max(smooth_dict[bar]) for i in smooth_dict[bar]]
        plt.plot(bar_ptime[bar][8:], smooth_dict[bar], linestyle=s, label=bar, color=colors[bar], alpha=a)
    plt.xlabel('Pseudotime')
    plt.ylabel('Avg of Seq'+str(seqlen_mod))
    plt.title('Model Importance of '+var)
    plt.legend(loc='best', fontsize=6)
    plt.ylim(bottom=0)
    plt.xlim(timezoom)
    plt.show()
    
    # only plot resistance
    for bar in smooth_dict.keys():
        if bar in resbars:
            s = 'dotted'
            plt.plot(bar_ptime[bar][8:], smooth_dict[bar], linestyle=s, label=bar, color=colors[bar])
    plt.xlabel('Pseudotime')
    plt.ylabel('Avg of Seq'+str(seqlen_mod))
    plt.title('Model Importance of '+var)
    plt.legend(loc='best', fontsize=6)
    plt.ylim(bottom=0)
    plt.xlim(timezoom)
    plt.show()
    
    # normalized to other features
    #plotFeatNormModel(barcodes, var, seqlen_mod, resbars, master, bar_ptime)
    
    # plot path scores from jdata with plotAvgFeaturePerBar
    plotAvgFeaturePerBar(var, jdata, barcodes, seqlen_act, resbars, timezoom)


# In[14]:


path = '/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/barmodels_combo_plsrres0.5_HALL3/'
master = masterBardict(new_barcodes, path)
#master[new_barcodes[0]]


# In[15]:


# top 20 paths
for p in countpathsl[:50]:
    print(p)


# In[25]:


path = '/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/barmodels_pcares0.25/'
master = masterBardict(new_barcodes, path)


# In[16]:


resbars9545 = ['GFPBC_libB_90850', 'GFPBC_libB_10678', 'GFPBC_libB_92196', 'GFPBC_libB_10737']
resbarspalbo = ['GFPBC_libB_90850', 'GFPBC_libB_92196']
resbarscombo = ['GFPBC_libB_37798', 'GFPBC_libB_92196']

ax = ridgePlot('')
for p in countpathsl[:30]:
    plotModelScoresPerBar(p[0], new_barcodes, data, 20, 50, master, resbarscombo, [0,1])


# In[17]:


def plotAsOne(var, barcodes, jdata, seqlen_mod, seqlen_act, master, resbars, timezoom, path):
    
    c=['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'yellow', 'black', 'yellow']
    colors = {barcodes[i]:c[i] for i in range(len(barcodes))}
    
    # get bardict for each barcode
    path_bardict = {}
    bar_ptime = {}
    for bar in barcodes:
        # lookup and store the error scores for path/gene
        path_bardict[bar] = master[bar][var]
        
        # get the ptime for cells in barcode from jdata
        df = pd.DataFrame(jdata[bar])
        bar_ptime[bar] = df['Pseudotime'].tolist()
    
    # smooth scores along ptime after normalized
    ##### mean of seqlen
    smooth_dict = {}
    for bar in path_bardict.keys():
        featbar = path_bardict[bar]
        seqnormscore = []
        for i in range(seqlen_mod, len(featbar)+seqlen_mod):
            seqnormscore.append(np.mean(featbar[i-seqlen_mod:i]))
        smooth_dict[bar] = seqnormscore
    
    
    # plot model data as max-norm plot
    fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(7.5, 5.5),
                        layout="constrained")
    
    
    for bar in smooth_dict.keys():
        if bar in resbars:
            s = 'dotted'
            a = 1
        else:
            s = 'solid'
            a = 0.25
        #featbar = [i / max(smooth_dict[bar]) for i in smooth_dict[bar]]
        axs[0].plot(bar_ptime[bar][8:], smooth_dict[bar], linestyle=s, label=bar, color=colors[bar], alpha=a)
    axs[0].set_xlabel('Pseudotime')
    axs[0].set_ylabel('Avg of Seq'+str(seqlen_mod))
    axs[0].set_title('Model Importance of '+var)
    axs[0].legend(loc='best', fontsize=6)
    axs[0].set_ylim(bottom=0)
    axs[0].set_xlim(timezoom)
    
    
    # plot path scores from jdata with plotAvgFeaturePerBar
    #plotAvgFeaturePerBar(var, jdata, barcodes, seqlen_act, resbars, timezoom)
    
    for b in barcodes:
        df = pd.DataFrame(jdata[b])
        featbar = df[var].tolist() ######## not normalized
        # normalize feature or not
        #featbar = [i / max(featbar) for i in featbar]
        ptime = df['Pseudotime'].tolist()
        ##### mean of 8-20 sequence
        seqnormscore = []
        for i in range(seqlen_act, len(featbar)+seqlen_act):
            seqnormscore.append(np.mean(featbar[i-seqlen_act:i]))
        if b in resbars:
            s = 'dotted'
            a = 1
        else:
            s = 'solid'
            a = 0.25
        barname = b.split('_')[-1]
        axs[1].plot(ptime, seqnormscore, linestyle=s, label=barname, color=colors[b], alpha=a)
        axs[1].set_xlabel('Pseudotime')
        axs[1].set_ylabel('Max-Normalized Avg of Seq'+str(seqlen_act))
        axs[1].set_title('Barcodes Activity of '+var)
        axs[1].set_xlim(timezoom)
        axs[1].legend(loc='best', fontsize=6)
        
    if path != '': 
        plt.savefig(path+var+'.png')
    plt.show()


# In[99]:


#plotAsOne('Score.HALLMARK_G2M_CHECKPOINT.1', new_barcodes, data, 20, 50, master, resbarspalbo, [0,1])


# In[20]:


path = '/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/barmodels_combo_plsrres0.5_HALL3/HALL/'

x = ridgePlot(path)
for p in countpathsl:
    plotAsOne(p[0], new_barcodes, data, 20, 50, master, resbarscombo, [0,1], path)


# In[54]:


def activity_vs_surv(var, barcodes, jdata, seqlen, master):
    
    # get bardict for each barcode
    path_bardict = {}
    bar_ptime = {}
    bar_surv = {}
    for bar in barcodes:
        # lookup and store the error scores for path/gene
        path_bardict[bar] = master[bar][var]
        
        # get the ptime for cells in barcode from jdata
        df = pd.DataFrame(jdata[bar])
        bar_ptime[bar] = df['Pseudotime'].tolist()
        bar_surv[bar] = df['Survival'].tolist()
    
    # smooth error and activity along ptime after normalized
    ##### mean of seqlen
    resbars = ['GFPBC_libB_90850', 'GFPBC_libB_10678', 'GFPBC_libB_92196', 'GFPBC_libB_10737']
    
    smooth_dict = {}
    xs = []
    ys = []
    for bar in path_bardict.keys():
        df_activity = pd.DataFrame(jdata[bar])
        featbar_activity = df_activity[var].tolist() ######## not normalized
        # to normalize
        #featbar_activity = [i / max(featbar_activity) for i in featbar_activity]
        ptime = df_activity['Pseudotime'].tolist()
        featbar_mod = path_bardict[bar]
        seqnormscore_mod = []
        seqnormscore_act = []
        #for i in range(seqlen, len(featbar_mod)+seqlen):
        #    seqnormscore_mod.append(np.mean(featbar_mod[i-seqlen:i]))
        #smooth_dict[bar] = seqnormscore_mod
        featbar_mod = [i / max(featbar_mod) if i > 0 else 0 for i in featbar_mod ]
        for i in range(seqlen, len(featbar_activity)+seqlen):
            seqnormscore_act.append(np.mean(featbar_activity[i-seqlen:i]))
        if bar in resbars:
            plt.scatter(bar_surv[bar], featbar_activity, color='black', alpha=featbar_mod, marker='^')
        if bar not in resbars:
            plt.scatter(bar_surv[bar], featbar_activity, color='black', alpha=featbar_mod, marker='o')
    plt.xlabel('Survival')
    plt.ylabel(var)
    plt.show()


# In[30]:





# In[14]:


os.chdir('/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/')

jfile1 = 'json_files/tran_9545_plsrres0.3_HALL3.json'
jfile2 = 'json_files/tran_palbo_plsrres0.3_HALL3.json'
jfile3 = 'json_files/tran_combo_plsrres0.5_HALL3.json'

import json
f1 = open(jfile1)
data1 = json.load(f1)

f2 = open(jfile2)
data2 = json.load(f2)

f3 = open(jfile3)
data3 = json.load(f3)


# In[15]:


def getTimeAndScore(feat, j, bars, seqlen):
    
    ptimedict = {}
    scoredict = {}
    
    for b in bars:
        df = pd.DataFrame(j[b])
        featbar = df[feat].tolist() ######## not normalized
        # normalize feature or not
        #featbar = [i / max(featbar) for i in featbar]
        ptime = df['Pseudotime'].tolist()
        ##### mean of 8-20 sequence
        seqnormscore = []
        for i in range(seqlen, len(featbar)+seqlen):
            seqnormscore.append(np.mean(featbar[i-seqlen:i]))
            
        ptimedict[b] = ptime
        scoredict[b] = seqnormscore
            
    return ptimedict, scoredict

def plotJdata(ptimedict, scoredict, treatment, colors, s, a):
    
    for b in ptimedict.keys():
        barname = b.split('_')[-1] + treatment
        plt.plot(ptimedict[b], scoredict[b], linestyle=s, label=barname, color=colors[b], alpha=a)
        
        

def plotComboScores(feat, jdata1, jdata2, jdata3, resbars1, resbars2, resbars3, seqlen, timezoom, treats):
    
    #c=['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'black', 'yellow', 'gray']
    c=['blue', 'orange', 'blue', 'orange', 'blue', 'orange']
    colors1 = {resbars1[i]:c[i] for i in range(len(resbars1))}
    colors2 = {resbars2[i]:c[i+len(resbars1)] for i in range(len(resbars2))}
    colors3 = {resbars3[i]:c[i+len(resbars1)+len(resbars2)] for i in range(len(resbars3))}
    
    ptimedict1, scoredict1 = getTimeAndScore(feat, jdata1, resbars1, seqlen)
    ptimedict2, scoredict2 = getTimeAndScore(feat, jdata2, resbars2, seqlen)
    ptimedict3, scoredict3 = getTimeAndScore(feat, jdata3, resbars3, seqlen)
    
    plotJdata(ptimedict1, scoredict1, treats[0], colors1, 'dotted', 0.5)
    plotJdata(ptimedict2, scoredict2, treats[1], colors2, 'dashed', 0.5)
    plotJdata(ptimedict3, scoredict3, treats[2], colors3, 'solid', 1)
    
    plt.xlabel('Pseudotime')
    plt.ylabel('Max-Normalized Avg of Seq'+str(seqlen))
    plt.title('Barcodes Activity of '+feat)
    plt.xlim(timezoom)
    plt.legend(loc='best', fontsize=6)
    plt.show()


# In[35]:


resbars9545 = ['GFPBC_libB_10678', 'GFPBC_libB_92196']
resbarspalbo = ['GFPBC_libB_10678', 'GFPBC_libB_92196']
resbarscombo = ['GFPBC_libB_10678', 'GFPBC_libB_92196']

plotComboScores('Score.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', data1, data2, data3, resbars9545, resbarspalbo, resbarscombo, 50, [0,1], ['9545', 'palbo', 'combo'])
plotComboScores('Score.HALLMARK_E2F_TARGETS.1', data1, data2, data3, resbars9545, resbarspalbo, resbarscombo, 50, [0,1], ['9545', 'palbo', 'combo'])
plotComboScores('Score.HALLMARK_ESTROGEN_RESPONSE_LATE', data1, data2, data3, resbars9545, resbarspalbo, resbarscombo, 50, [0,1], ['9545', 'palbo', 'combo'])
plotComboScores('Score.HALLMARK_G2M_CHECKPOINT.1', data1, data2, data3, resbars9545, resbarspalbo, resbarscombo, 50, [0,1], ['9545', 'palbo', 'combo'])
plotComboScores('Score.TCGA_LuBrCa_pAKT_dn', data1, data2, data3, resbars9545, resbarspalbo, resbarscombo, 50, [0,1], ['9545', 'palbo', 'combo'])
for p in countpathsl[:40]:
    plotComboScores(p[0], data1, data2, data3, resbars9545, resbarspalbo, resbarscombo, 50, [0,1], ['9545', 'palbo', 'combo'])


# In[22]:


resbars9545 = ['GFPBC_libB_90850', 'GFPBC_libB_10678']
resbarspalbo = ['GFPBC_libB_90850', 'GFPBC_libB_10678']
resbarscombo = ['GFPBC_libB_90850', 'GFPBC_libB_10678']

plotComboScores('Score.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', data1, data2, data3, resbars9545, resbarspalbo, resbarscombo, 50, [0,1], ['9545', 'palbo', 'combo'])
plotComboScores('Score.HALLMARK_E2F_TARGETS.1', data1, data2, data3, resbars9545, resbarspalbo, resbarscombo, 50, [0,1], ['9545', 'palbo', 'combo'])
plotComboScores('Score.HALLMARK_ESTROGEN_RESPONSE_LATE', data1, data2, data3, resbars9545, resbarspalbo, resbarscombo, 50, [0,1], ['9545', 'palbo', 'combo'])
plotComboScores('Score.HALLMARK_G2M_CHECKPOINT.1', data1, data2, data3, resbars9545, resbarspalbo, resbarscombo, 50, [0,1], ['9545', 'palbo', 'combo'])
plotComboScores('Score.TCGA_LuBrCa_pAKT_dn', data1, data2, data3, resbars9545, resbarspalbo, resbarscombo, 50, [0,1], ['9545', 'palbo', 'combo'])


# In[32]:


os.chdir('/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/')
jfile = 'json_files/tran_9545_plsrres0.3_HALL3.json'

import json
f = open(jfile)
data = json.load(f)

barcodes = [i for i in data.keys()]
new_barcodes = []
for i in barcodes:
    df = pd.DataFrame(data[i])
    if df.shape[0] > 200: # 50 or 60
        new_barcodes.append(i)

path = '/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/barmodels_9545_plsrres0.3_HALL3/'
master = masterBardict(new_barcodes, path)


var = 'barmodels_9545_plsrres0.3_HALL3'
os.chdir('/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/'+var)
feat = 'Score.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'
#barcode = 'GFPBC_libB_92196'
seqlen = 50
seqlen_mod = 20


# In[33]:


import matplotlib.pyplot as plt

barcode = 'GFPBC_libB_92196'
var = 'barmodels_9545_plsrres0.3_HALL3'
os.chdir('/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/'+var)

df = pd.read_csv('transformer_full'+barcode+'_feature_importance_full.csv')
basedf = df[df['feature'] == 'BASELINE']
realp = eval(basedf['prediction'].tolist()[0])
realt = eval(basedf['true'].tolist()[0])
featerr = {}
smalldf = df[df['feature'] == feat]


pdf = pd.DataFrame()
for col in smalldf.columns.tolist():
    if 'prediction_' in col:
        pdf[col] = eval(smalldf[col].tolist()[0])

pdf['Real_Prediction'] = realp
pdf['Real_True'] = realt

predictions = []
for each in pdf.index.tolist():
    predictions.append(np.mean(pdf.iloc[each].tolist()[:30]))

seqnormscore = []
for i in range(seqlen, len(predictions)+seqlen):
    seqnormscore.append(np.mean(predictions[i-seqlen:i]))

os.chdir('/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/')
import json
f = open(jfile)
data = json.load(f)
df = pd.DataFrame(data[barcode])
ptime = df['Pseudotime'].tolist()


plt.figure(figsize=(20, 5), dpi=80)
plt.plot(ptime[8:], realt, color='black', alpha=0.75, label='True Survival')
plt.plot(ptime[8:], seqnormscore, color='red', alpha=0.5, label='Permuted '+feat)
#plt.plot(range(len(realp)), realp, color='blue')
plt.ylabel('Normalized Clonal Survival')
plt.xlabel('Pseudotime')
plt.title('Survival and '+feat+' Feature Permutation for '+barcode)
plt.legend(loc='best')
plt.show()



# get bardict for each barcode
path_bardict = {}
bar_ptime = {}

# lookup and store the error scores for path/gene
path_bardict[barcode] = master[barcode][feat]

# smooth scores along ptime after normalized
##### mean of seqlen
smooth_dict = {}
featbar = path_bardict[barcode]
seqnormscore = []
for i in range(seqlen_mod, len(featbar)+seqlen_mod):
    seqnormscore.append(np.mean(featbar[i-seqlen_mod:i]))
smooth_dict[barcode] = seqnormscore

seqnormscore2 = []
for i in range(150, len(featbar)+150):
    seqnormscore2.append(np.mean(smooth_dict[barcode][i-150:i]))


plt.figure(figsize=(20, 5), dpi=80)
for bar in smooth_dict.keys():
    plt.plot(ptime[8:], smooth_dict[barcode], linestyle='solid', alpha=0.25, label=barcode+' avg of 20', color='blue')
    plt.plot(ptime[8:], seqnormscore2, linestyle='solid', alpha=1, label=barcode+' avg of 150', color='blue')
plt.xlabel('Pseudotime')
plt.ylabel('Smoothed Model Importance')
plt.title('Model Importance of '+barcode+' '+feat)
plt.legend(loc='best', fontsize=6)
plt.ylim(bottom=0)
#plt.xlim(timezoom)
plt.show()


# In[ ]:




