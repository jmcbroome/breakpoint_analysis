#!/usr/bin/env python3

#import
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random
import statistics as st
from chromatin_script import chromatin, pixelate, tandem_dup, getchrom, sumchrom, makedist, chrombinary, GFFRead
#define functions/classes
chromlen = {'2L':23515712, '2R':25288936, '3L':25288936, '3R':32081331, 'X':23544271, '4':1830000} #values as of DM6 assembly, chr 2/3/X only for first testing, chr4 an estimate

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    #add args
    parser.add_argument('-f', '--file', help = 'path to input inversion data.', default = '../common_inv_all.txt')
    parser.add_argument('-c', '--chromatin', help = 'path to input chromatin data.', default = '../droso_chromatin_r6.gff3')
    parser.add_argument('-x', '--fixed', help = 'path to fixed inversion data.', default = '../fixed_inv_2.txt')
    args = parser.parse_args()
    return args

def build_expected(crdata, pixels):
    def makefake(crdata, pixels, factor = 10):
        fakecr = []
        for entry in np.array(crdata):
            for f in range(0,factor):
                nent = [e for e in entry]
                nent[2] = random.choice(range(0,chromlen[entry[0]] - abs(int(entry[3])-int(entry[2]))))
                nent[3] = nent[2] + abs(int(entry[3])-int(entry[2]))
                fakecr.append(nent)
        fcr = pd.DataFrame(fakecr)
        fcr.columns = crdata.columns
        chrodist= makedist(np.array(fcr), pixels)
        fcr.loc[:,range(1,10)] = chrodist
        trinvec = chrombinary(np.array(fcr.loc[:,range(1,10)]))
        rerep = {k:0 for k in set(trinvec)}
        for v in trinvec:
            rerep[v] += 1
        trinvec = pd.DataFrame(trinvec)
        trinvec.columns = ["Activity"]
        fcr['Activity'] = trinvec['Activity']
        return fcr

    perm_exp = {i:[] for i in range(0,9)}
    for p in range(0,1000):
        fcr = makefake(crdata, pixels)
        farrays = {'freq':[], 'act':[], 'count':[]}
        for ftype in ['Rare','Common','Fixed']:
            bpdat = fcr.loc[fcr['Freq'] == ftype]
            for atype in ['Active','Mixed','Inactive']:
                bpdat2 = bpdat.loc[bpdat['Activity'] == atype]
                farrays['freq'].append(ftype)
                farrays['act'].append(atype)
                farrays['count'].append(len(bpdat2))
        #divide all numbers by the factor 10 to get the expected values of each category
        exp = []
        for e in farrays['count']:
            exp.append(e/10)

        for i,e in enumerate(exp):
            perm_exp[i].append(e)
    #use means to build dif, calculate standard deviations
    #calcualate differents and add to arrays object, then graph the differences
    edif = []
    stds = []
    for i,v in perm_exp.items():
        mean_exp = st.mean(v)
        edif.append(mean_exp)
        stds.append(st.stdev(v))

    #convert to differences from the mean?
    perm_dif = {}
    means = []
    for i,v in perm_exp.items():
        mean = st.mean(v)
        means.append(mean)
        nv = []
        for p in v:
            nv.append(p-mean)
        perm_dif[i] = nv
    return perm_dif, means
    # #apply data labeling
    # cats = {0:['Rare','Active'],
    #         1:['Rare','Mixed'],
    #         2:['Rare','Inactive'],
    #         3:['Common','Active'],
    #         4:['Common','Mixed'],
    #         5:['Common','Inactive'],
    #         6:['Fixed','Active'],
    #         7:['Fixed','Mixed'],
    #         8:['Fixed','Inactive']
    #     }
    # freqv=[]
    # actv=[]
    # countv=[]
    # for i,v in perm_dif.items():
    #     for p in v:
    #         freqv.append(cats[i][0])
    #         actv.append(cats[i][1])
    #         countv.append(p)
    # return freqv, countv, actv

def draw_box(xval, ydist, panel, width, color = 'black'):
    order = sorted(ydist)
    fifth = np.percentile(order, 5)
    twentyfifth = np.percentile(order, 25)
    half = np.percentile(order, 50)
    seventyfifth = np.percentile(order, 75)
    ninetyfifth = np.percentile(order, 95)
    panel.plot([xval-width/2,xval+width/2],[fifth,fifth], color = color, linewidth = 1)
    panel.plot([xval,xval],[fifth, twentyfifth], color = color, linewidth = 1)
    panel.plot([xval-width,xval+width],[twentyfifth,twentyfifth], color = color, linewidth = 1)
    panel.plot([xval-width,xval-width],[twentyfifth,half], color = color, linewidth = 1)
    panel.plot([xval+width,xval+width],[twentyfifth,half], color = color, linewidth = 1)
    panel.plot([xval-width,xval+width],[half,half], color = color, linewidth = 1)
    panel.plot([xval-width,xval-width],[half, seventyfifth], color = color, linewidth = 1)
    panel.plot([xval+width,xval+width],[half, seventyfifth], color = color, linewidth = 1)
    panel.plot([xval-width,xval+width],[seventyfifth,seventyfifth], color = color, linewidth = 1)
    panel.plot([xval,xval],[seventyfifth, ninetyfifth], color = color, linewidth = 1)
    panel.plot([xval-width/2,xval+width/2],[ninetyfifth,ninetyfifth], color = color, linewidth = 1)
    return panel

def main():
    args = argparser()
    
    #read in the data
    with open(args.file, 'r') as crinf:
    
        #read in the .tsv file with inversion coordinates and frequencies
        common_rares = []
        for entry in crinf.readlines()[1:]:
            spent = entry.strip().split()
            try:
                if len(spent) >= 9:
                    invid = spent[0]
                    chrom = invid[3:5]
                    rare = spent[1]
                    if chrom[1] == ')':
                        chrom = chrom[0]
                    if spent[3][0] == 'P' or spent[3][0] == 'C':
                        common_rares.append([chrom, invid, spent[4], spent[5], rare])
                    else:
                        common_rares.append([chrom, invid, spent[3], spent[4], rare])
                else:
                    common_rares.append([chrom, invid, spent[0], spent[1], rare])
            except IndexError:
                continue
        #define a dataframe that will contain the break data for the rest of the analysis
    crdata = pd.DataFrame(common_rares[:-3])
    crdata.columns = ['Chro','Label','Forward','Reverse','Freq']
    dups = tandem_dup(crdata)
    dups = pd.DataFrame(dups)
    dups.columns = ["TanDup"]
    crdata = pd.concat([crdata, dups], 1)
    fixed_singles = []
    with open(args.fixed, 'r') as fixinv:
        for entry in fixinv.readlines():
            spent = entry.strip().split()
            #need a unique identifier for each one.
            if len(spent) == 9:
                invid = spent[0]
                fixed_singles.append([spent[1], invid, spent[6], spent[7], spent[8]])
            else:
                fixed_singles.append([spent[0], invid, spent[5], spent[6], spent[7]])
    fixed_singles = pd.DataFrame(fixed_singles)
    fixed_singles.columns = ['Chro','Label','Forward','Reverse','TanDup']
    fixed_singles['Freq'] = 'Fixed'
    fixed_singles = fixed_singles.drop(19)
    crdata = pd.concat([crdata, fixed_singles], ignore_index = True, sort = True)
    crdata = crdata[['Chro','Label','Forward','Reverse','Freq','TanDup']]
    
    gff = GFFRead(args.chromatin)[1]
    nines = []
    for entry in gff:
        if int(entry[4][3]) == 9:
            nines.append([entry[0][:2]] + entry[1:4] + ['9_state'])
    #analysis proceeds with only nine-state data from this point; uses nine-state to determine active/inactive/mixed status of region.    
    test = chromatin(nines)
    pixels = pixelate(test.chrom, 10000)
    chrodist= makedist(np.array(crdata), pixels)
    crdata = pd.concat([crdata, chrodist], 1)

    trinvec = chrombinary(np.array(crdata.loc[:,range(1,10)]))
    # print(trinvec)
    trinvec = pd.DataFrame(trinvec)
    trinvec.columns = ["Activity"]
    crdata = pd.concat([crdata, trinvec], 1)
    if args.verbose:
        print("Inversion and Chromatin Data read")

    pdifs, means = build_expected(crdata, pixels)
    arrays = {'freq':[], 'act':[], 'count':[]}
    for ftype in ['Rare','Common','Fixed']:
        bpdat = crdata.loc[crdata['Freq'] == ftype]
        for atype in ['Active','Mixed','Inactive']:
            bpdat2 = bpdat.loc[bpdat['Activity'] == atype]
            arrays['freq'].append(ftype)
            arrays['act'].append(atype)
            arrays['count'].append(len(bpdat2))
    nvs = []
    # typeconv = {'Rare':0,'Common':1,"Fixed":2}
    for i in range(len(arrays['count'])):
        nv = arrays['count'][i] - means[i]
        nvs.append(nv)
    if args.verbose:
        print("Permutations Complete")
    #pdifs is indexed with values 0-8
    #0=Rare/Active, 1=Rare/Mixed, 2=Rare/Inactive... etc until 8=Fixed/Inactive (rare-common-fixed/active-mixed-inactive)

    #set up figure structure
    plt.figure(figsize=[10,8])
    panel1 = plt.axes([1/5,1/4,3/5,1/2])
    plt.title('Chromatin State of Inversion Breakpoints')
    panel1.set_ylabel('Difference')
    panel1.set_ylim(-6,10)
    panel1.set_yticks(list(range(-6,10,2)))
    panel1.set_xlim(-.5,2.5)
    panel1.set_xticks([0,1,2])
    panel1.set_xticklabels(["Rare","Common","Fixed"])
    colors = {0:'red',1:'purple',2:'blue',
        3:'red',4:'purple',5:'blue',
        6:'red',7:'purple',8:'blue'}
    for xval in [0,1,2]:
        for i in [0,1,2]:
            pdi = xval*3+i
            panel1 = draw_box(xval,pdifs[pdi],panel1,2/7,colors[pdi])
            panel1.scatter(x=[xval],y=[nvs[pdi]],s=20,color=colors[pdi])
    
    # plt.autoscale()
    
    # plt.show()
    plt.savefig('figure_s3.png',dpi=200)

if __name__ == "__main__":
    main()