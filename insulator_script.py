#!/usr/bin/env python3

#import
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statistics as st
from scipy.stats import mannwhitneyu, fisher_exact
import math
import random
#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-f', '--file', help = 'Path to the file with common and rare breakpoint information', default = 'common_inv_all.txt')
    parser.add_argument('-x', '--fixed', help = 'Path to the file containing fixed breakpoints (different format from the other file', default = 'fixed_inv_2.txt')
    parser.add_argument('-i', '--insulator', help = 'Path to the insulator annotation file', default = 'public_datasets/GSE26905_Dm_Insulator_Classes_r6.GFF3')
    args = parser.parse_args()
    return args

def tandem_dup(crdata):
    dups = []
    for entry in np.array(crdata):
        dups.append(abs(int(entry[2])-int(entry[3])))
    return dups

def insulator_search(data, chrom, start, end): #duplicated for reference.
    #just an iterator because the file isn't that big anyways.
    sites = []
    for entry in data[chrom]:
        if start < entry < end:
            sites.append(entry)
    return sites

def insulator_track(classx, bpdata):
    dvec = []
    for spent in np.array(bpdata):
        minis = 100000000000000
        minie = 100000000000000
        for element in classx[spent[0]]:
            if abs(int(element) - int(spent[2])) < abs(minis):
                minis = abs(int(element) - int(spent[2]))
                if int(spent[2]) <= int(element) < int(spent[3]):
                    sign = -1
                else:
                    sign = 1
            if abs(int(element) - int(spent[3])) < abs(minie):
                minie = abs(int(element)- int(spent[3]))
                if int(spent[2]) <= int(element) < int(spent[3]):
                    sign = -1
                else:
                    sign = 1
        if abs(minis) < abs(minie):
            dvec.append(abs(minis*sign)) #reverse sign of forward breaks, because I care about outside duplicate as positive and inside as negative
        else:
            dvec.append(abs(minie*sign)) 
    return dvec
    
def insulator_permuter(classx, bpdata, pernum = 100):
    pervec = []
    for entry in np.array(bpdata):
        sim = []
        for p in range(0,pernum):
            nent = [e for e in entry]
            nent[2] = random.choice(range(0,chromlen[entry[0]] - abs(int(entry[3])-int(entry[2]))))
            nent[3] = nent[2] + abs(int(entry[3])-int(entry[2]))
            sim.append(insulator_track(classx, [nent])[0])
#             sim.append(abs(insulator_track(classx, [nent])[0])) #ignore sign of distance for now.
        #then do a statistical check whether this particular inversion breakpoint is nearer or farther from an insulator than expected at random.
        #sim is the distribution of expected distances with length pernum
        #since its supposed to be closer, lets check whether its shorter than the 5% place in the sorted sim list.
        rd = insulator_track(classx, [entry])[0]
        pervec.append(sim)
    
    return pervec

def main():
    args = argparser()
    #insert code    
    common_rares = []
    with open(args.file, 'r') as cinv:
        for entry in cinv.readlines()[1:]:
            spent = entry.strip().split()
            try:
                if len(spent) >= 9:
                    invid = spent[0]
                    chrom = invid[3:5]
                    rare = spent[1]
    #                 if chrom[0] == '7':
    #                     print(spent)
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
    crdata = pd.concat([crdata, fixed_singles], ignore_index = True)
    crdata = crdata[['Chro','Label','Forward','Reverse','Freq','TanDup']]
    #read in insulators
    with open(args.insulator, 'r') as insulators:
        #data structure is two dicts, key of chromosome, value of a set of insulator element points
        #treating as point bases because all binding sites are 10 bases long
        chromlen = {'2L':23515712, '2R':25288936, '3L':25288936, '3R':32081331, 'X':23544271, '4':1830000} #values as of DM6 assembly, chr 2/3/X only for first testing, chr4 an estimate
        class1 = {k:[] for k in chromlen.keys()}
        class2 = {k:[] for k in chromlen.keys()}
        for entry in insulators.readlines():
            spent = entry.strip().split()
            if len(spent) > 8:
                try:
                    if spent[9] == 'I':
                        class1[spent[0][3:]].append(int(spent[3]))
                    elif spent[9] == "II":
                        class2[spent[0][3:]].append(int(spent[3]))
                except:
                    continue
    #class 1 insulators are the ones that exhibit blockage properties in the literature
    t = insulator_permuter(class1, crdata)
    td = pd.DataFrame(columns = ['InsDist_exp'])
    for i,it in enumerate(t):
        td.loc[i] = st.median(it)
    crdata = pd.concat([crdata,td],1)
    crdata = crdata.loc[:,~crdata.columns.duplicated()]
    realvec = [v for v in insulator_track(class1, crdata)]
    rv = pd.DataFrame(realvec)
    rv.columns = ['InsDist']
    crdata = pd.concat([crdata,rv],1)
    crdata = crdata.loc[:,~crdata.columns.duplicated()]

    c = [e for e in crdata.loc[crdata['Freq'] == 'Common']['InsDist'] if abs(e) < 100000]
    a = [e for e in crdata.loc[crdata['Freq'] == 'Rare']['InsDist'] if abs(e) < 100000]
    x = [e for e in crdata.loc[crdata['Freq'] == 'Fixed']['InsDist'] if abs(e) < 100000]
    # f = [e for e in crdata['InsDist_exp'].astype(int) if abs(e) < 100000]
    fperms = insulator_permuter(class1, crdata, 1000)
    #flatten fperms
    f = []
    for fs in fperms:
        for s in fs:
            f.append(s)
    f = [e for e in f if abs(e) < 100000]
    print('Fixed closer than Expected:',mannwhitneyu(x,f,alternative='less')[1])
    print('Coommon closer than Expected:',mannwhitneyu(c,f,alternative='less')[1])
    print('Rare closer than Expected:',mannwhitneyu(a,f,alternative='less')[1])
    print('Fixed closer than Common:',mannwhitneyu(c,x,alternative='greater')[1])
    print('Fixed closer than Rare:',mannwhitneyu(a,x,alternative='greater')[1])
    print('Common closer than Rare:',mannwhitneyu(c,a,alternative='less')[1])


if __name__ == "__main__":
    main()