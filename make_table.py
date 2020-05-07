#NOT A COMPLETE IMPLEMENTATION.
#import
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statistics as st
from scipy.stats import mannwhitneyu
import math
import random
#import from other scripts.
from genome_annote_script import inter_checker
from insulator_script import insulator_search, insulator_track
from tandup_script import tandem_dup
from chromatin_script import *

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-f', '--file', help = 'Path to the file with common and rare breakpoint information', default = '../common_inv_all.txt')
    parser.add_argument('-x', '--fixed', help = 'Path to the file containing fixed breakpoints (different format from the other file', default = '../fixed_inv_2.txt')
    parser.add_argument('-c', '--chromatin', help = 'Path to the chromatin file', default = '../public_datasets/droso_chromatin_r6.GFF3')
    parser.add_argument('-a', '--annote', help = 'Path to the annotation file.', default = '../public_datasets/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff')
    parser.add_argument('-p', '--polytene', help = 'Path to the polytene TAD file.', default = '../public_datasets/POLYTENE_TADS_R6.txt')
    parser.add_argument('-i', '--insulator', help = 'Path to the insulator annotation file', default = '../public_datasets/GSE26905_Dm_Insulator_Classes_r6.GFF3')
    parser.add_argument('-t', '--tads', help = 'Path to the non-polytene TAD file', default = '../public_datasets/sexton_tads_fulldata.txt')
    args = parser.parse_args()
    return args

chromlen = {'2L':23515712, '2R':25288936, '3L':28110227, '3R':32081331, 'X':23544271, '4':1830000} #values as of DM6 assembly, chr 2/3/X only for first testing, chr4 an estimate

def read_annotation(args):
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
    return annotation

def read_sexton(args, annotation):
    #additional tad data incorporated as a separate set
    sexton_tad = {k:[] for k in chromlen.keys()}
    with open(args.tads,'r') as ptadf:
        for entry in ptadf.readlines():
            if len(entry.strip()) > 0:
                spent = entry.strip().split()
                if spent[0] == '"Table' or spent[0] == 'domain':
                    continue
                sexton_tad[spent[1]].append((int(spent[2]), int(spent[3]), spent[0] + "_" + spent[4]))
    annotation['sexton_tad'] = sexton_tad
    return annotation

def main():
    args = argparser()
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
    #define the dataframe.
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
    fixed_singles = fixed_singles.drop(19) #this particular break is bad data
    crdata = pd.concat([crdata, fixed_singles], ignore_index = True)
    crdata = crdata[['Chro','Label','Forward','Reverse','Freq']]
    crdata['TanDup'] = [abs(int(crdata['Forward'][i])-int(crdata['Reverse'][i])) for i in range(len(crdata['Forward']))]
    #add chromatin data
    gff = GFFRead(args.chromatin)[1]
    nines = []
    for entry in gff:
        if int(entry[4][3]) == 9:
            nines.append([entry[0][:2]] + entry[1:4] + ['9_state'])
    test = chromatin(nines)
    pixels = pixelate(test.chrom, 10000)
    chrodist= makedist(np.array(crdata), pixels)
    crdata = pd.concat([crdata, chrodist], 1)
    trinvec = chrombinary(np.array(crdata.loc[:,range(1,10)]))
    crdata['Activity'] = trinvec
    #add insulator data
    with open(args.insulator, 'r') as insulators:
        #data structure is two dicts, key of chromosome, value of a set of insulator element points
        #treating as point bases because all binding sites are 10 bases long
        chromlen = {'2L':23515712, '2R':25288936, '3L':28110227, '3R':32081331, 'X':23544271, '4':1830000} #values as of DM6 assembly, chr 2/3/X only for first testing, chr4 an estimate
        class1 = {k:[] for k in chromlen.keys()}
        for entry in insulators.readlines():
            spent = entry.strip().split()
            if len(spent) > 8:
                try:
                    if spent[9] == 'I':
                        class1[spent[0][3:]].append(int(spent[3]))
                except:
                    continue
    crdata['InsDist'] = [v for v in insulator_track(class1, crdata)]
    #add annotation interruption data
    annotation = read_annotation(args)
    annotation = read_sexton(args, annotation) #additional annotation data incorporated.
    int_vecs = {e:[] for e in ['gene','mRNA','polytad']}
    for etype, chroms in annotation.items():
        if etype in int_vecs.keys():
            for i, entry in crdata.iterrows():
                state = inter_checker((entry['Chro'],entry['Forward'],entry['Reverse']),chroms)
                int_vecs[etype].append(state)
    for e, vec in int_vecs.items():
        crdata[e] = vec
    #more complex sexton et al tad data
    chroms = annotation['sexton_tad']
    sexton_int = []
    sexton_lab = []
    sexton_state = []
    for i, entry in crdata.iterrows():
        ints, edat = inter_checker((entry['Chro'],entry['Forward'],entry['Reverse']),chroms, etype = True)
        sexton_int.append(ints)
        label, state = edat.split("_")
        sexton_lab.append(label)
        sexton_state.append(state)
    crdata['sexton_tad'] = sexton_int
    crdata['sexton_lab'] = sexton_lab
    crdata['sexton_state'] = sexton_state
    crdata.to_csv('full_crdata.csv')

if __name__ == "__main__":
    main()    