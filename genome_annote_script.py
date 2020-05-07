#!/usr/bin/env python3

#import
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statistics as st
from scipy.stats import mannwhitneyu, fisher_exact, percentileofscore
import math
import random

chromlen = {'2L':23515712, '2R':25288936, '3L':28110227, '3R':32081331, 'X':23544271, '4':1830000} #values as of DM6 assembly, chr 2/3/X only for first testing, chr4 an estimate

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-f', '--file', help = 'Path to the file with common and rare breakpoint information', default = '../common_inv_all_v2.tsv')
    parser.add_argument('-x', '--fixed', help = 'Path to the file containing fixed breakpoints (different format from the other file', default = '../fixed_inv_2.txt')
    parser.add_argument('-a', '--annote', help = 'Path to the annotation file', default = '../public_datasets/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff')
    parser.add_argument('-p', '--polytene', help = 'Path to the polytene TAD file.', default = '../public_datasets/POLYTENE_TADS_R6.txt')
    parser.add_argument('-t', '--tads', help = 'Path to the non-polytene TAD file.', default = '../public_datasets/sexton_tads_r6.txt')
    args = parser.parse_args()
    return args
    
def tandem_dup(crdata):
    dups = []
    for entry in np.array(crdata):
        dups.append(abs(int(entry[2])-int(entry[3])))
    return dups

def inter_checker(breakpoint, annote, etype = False):
    #takes a single breakpoint coordinate as a tuple of chro and two integers and checks to see what it interrupts in annote
    #if etype == True, returns the fourth string value of the tuple. 
    for element in annote[breakpoint[0]]:
        if element[0] < int(breakpoint[1]) < element[1] and element[0] < int(breakpoint[2]) < element[1]:
            if etype:
                return ("Int", element[2])
            else:  
                return "Int"
    if etype:
        return ('Non', 'N/A_N/A')
    else:
        return "Non"

def inter_permuter(crdata, annotations, pnum = 500):
    einterrat = {'Non':0, 'Int':0}
    subsamples = []
    for p in range(0, pnum):
        sub = {k:0 for k in einterrat.keys()}
        for i, entry in crdata.iterrows():
            rand = random.choice(range(0,chromlen[entry['Chro']]-int(entry['TanDup'])))
            fbp = (entry['Chro'], rand, rand + int(entry['TanDup']))
            fint = inter_checker(fbp, annotations)
            einterrat[fint] += 1
            sub[fint] += 1
        subsamples.append([sub[e] for e in ['Non', 'Int']])
    rinterrat = {'Non':0, 'Int':0}
    for i, entry in crdata.iterrows():
        rbp = (entry['Chro'], entry['Forward'], entry['Reverse'])
        rint = inter_checker(rbp, annotations)
        rinterrat[rint] += 1
    pval = fisher_exact([list(einterrat.values()),list(rinterrat.values())])[1]
    return einterrat, rinterrat, pval, subsamples

def main():
    args = argparser()
    
    annotation = {}
    genespan = {}
    cchro = 'X'
    cgene = 'none'
    ingene = False
    with open(args.annote,'r') as fullannote:
        for entry in fullannote.readlines():
    #         try:
            #this is a friggin complicated file.
            #rows I care about have id,RefSeq, element type (!), start, stop, some kind of code with -,+,- for 3 columns, then a super long ID I'd need to regex through.
            #some rows start with hashes, skip those
            #the chromosome everything is on depends on the location in the file. 
            spent = entry.strip().split()
            if spent[0][0] == "#":
                continue #skip hash rows.
            if spent[2] == 'region':
                #here's where parsing gets complicated.
                string = spent[8] #the super long nasty thing.
                #I'm going to want to look along string until I locate the entry 'chromosome='
                #then grab everything after that up to but not including a ';'
                c = []
                for i in range(0, len(string)-11):
                    if string[i:i+11] == 'chromosome=':
                        ci = i + 11
                        while string[ci] != ';':
                            c.append(string[ci])
                            ci += 1
                cchro = ''.join(c)
            else:
                annotation[spent[2]] = annotation.get(spent[2], {})
                annotation[spent[2]][cchro] = annotation[spent[2]].get(cchro, []) + [(int(spent[3]), int(spent[4]))]
                if spent[2] == 'exon':
                    ingene = True
                    genespan[cchro] = genespan.get(cchro, [])
                    string = spent[8]
                    for i in range(0,len(string) - 14):
                        if string[i:i+7] == 'GeneID:':
                            i += 7
                            if string[i:i+7] != cgene:
                                cgene = string[i:i+7]
                                genespan[cchro].append([int(spent[3])])
                else:
                    if ingene:
                        cont = False
                        for i in range(0,len(string) - 14):
                            if string[i:i+7] == 'GeneID:':
                                i += 7
                                if string[i:i+7] == cgene:
                                    cont = True
                        if not cont:
                            ingene = False
                            try:
                                # print(genespan[cchro][-1].append(int(spent[4])))
                                genespan[cchro][-1].append(int(spent[4]))
                            except KeyError:
                                continue
                                # print(cchro)
                                # print(genespan.keys())

    polytad = {k:[] for k in chromlen.keys()}
    with open(args.polytene,'r') as ptadf:
        for entry in ptadf.readlines():
            spent = entry.strip().split()
            polytad[spent[0]].append((int(spent[1]), int(spent[2])))
    annotation['polytad'] = polytad
    #additional tad data incorporated as a separate set
    sexton_tad = {k:[] for k in chromlen.keys()}
    with open(args.tads,'r') as ptadf:
        for entry in ptadf.readlines():
            spent = entry.strip().split()
            sexton_tad[spent[0]].append((int(spent[1]), int(spent[2])))
    annotation['sexton_tad'] = sexton_tad

    common_rares = []
    with open(args.file, 'r') as cinv:
        for entry in cinv.readlines()[1:]:
            spent = entry.strip().split()
            try:
                if (len(spent) >= 9 and spent[1] == 'Common') or len(spent) >= 10:
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
    fixed_singles = []
    with open(args.fixed, 'r') as fixinv:
        for entry in fixinv.readlines():
            spent = entry.strip().split()
            #need a unique identifier for each one.
            if len(spent) == 9:
                invid = spent[0]
                fixed_singles.append([spent[1], invid, spent[6], spent[7]])
            else:
                fixed_singles.append([spent[0], invid, spent[5], spent[6]])
    fixed_singles = pd.DataFrame(fixed_singles)
    fixed_singles.columns = ['Chro','Label','Forward','Reverse']
    fixed_singles['Freq'] = 'Fixed'
    fixed_singles = fixed_singles.drop(19)
    crdata = pd.concat([crdata, fixed_singles], ignore_index = True)
    crdata = crdata[['Chro','Label','Forward','Reverse','Freq']]
    crdata['TanDup'] = [abs(int(crdata['Forward'][i])-int(crdata['Reverse'][i])) for i in range(len(crdata['Forward']))]

    #drop the potential false positive rares.
    #names = ['In(2R)DL4', 'In(3L)DL5', 'In(3L)DL7', 'In(3L)DL9', 'In(3L)DL11', 'In(3R)DL12', 'In(3R)DL13']
    #crdata = crdata.loc[~crdata['Label'].isin(names)].reset_index(drop=True)
    #or drop all but high confidence rares
    #names = ['In(2L)DL3','In(2R)Y1a','In(2R)MAL_1*','In(2R)MAL_2*',"In(3L)DL10","In(3R)Gb"]
    #crdata = crdata.loc[(crdata['Label'].isin(names)) | (crdata['Freq'] != 'Rare')].reset_index(drop=True)    
    
    #run tests
    # gathered_subsamples = {f:{e:[] for e in annotation.keys()} for f in list(set(crdata['Freq']))}
    # gathered_rvals = {f:{e:[] for e in annotation.keys()} for f in list(set(crdata['Freq']))}
    #now do some comparatives for gene and mrna based stuff.
#now do some comparatives for gene and mrna based stuff.
    for etype in ['gene','polytad']:
    # for etype in ['polytad']:
        print("Analyzing {}".format(etype))
        tdat = {}
        for ftype in ['Fixed', 'Common', 'Rare', 'ComRare']:
            if ftype != 'ComRare':
                crdat1 = crdata.loc[crdata['Freq'] == ftype]
            else:
                crdat1 = crdata[crdata.Freq.isin(['Common','Rare'])]
    #         for ftype2 in ['Fixed','Common','Rare']:
    #             crdat2 = crdata.loc[crdata['Freq'] == ftype2]
    #             if ftype != ftype2:
            e1,r1,p1,ss1 = inter_permuter(crdat1,annotation[etype], pnum = 1000)
            tdat[ftype] = (e1, r1, p1, ss1) 
    #         e2,r2,p2,ss2 = inter_permuter(crdat2,annotation[etype], pnum = 1000, adjust = False)
        for k, vs in tdat.items():
            e1, r1, p1, ss1 = vs
            print("For {}, real value {}, percentile {}, where {} is 5th percentile of expectations and {} is 95th percentile of expectations".format(k, r1['Int'], percentileofscore([e[1] for e in ss1], r1["Int"]), np.percentile([e[1] for e in ss1], 5), np.percentile([e[1] for e in ss1],95)))
            fe_pv = fisher_exact([[r1['Non'], r1["Int"]],[np.mean([e[0] for e in ss1]), np.mean([e[1] for e in ss1])]])[1]
            print([list(r1.values()),[np.mean([e[0] for e in ss1]), np.mean([e[1] for e in ss1])]])
            print("Fisher's Exact Comparison pval: {}".format(fe_pv))
            print("##") 
if __name__ == "__main__":
    main()