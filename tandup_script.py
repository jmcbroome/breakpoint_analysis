#!/usr/bin/env python3

#import
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statistics as st
from scipy.stats import mannwhitneyu
import random
import math

#define functions/classes
def tandem_dup(crdata):
    '''
    Calculate tandem duplication sizes across the dataframe.
    '''
    dups = []
    for entry in np.array(crdata):
        dups.append(abs(int(entry[2])-int(entry[3])))
    return dups

#setting up what I will need to make the swarmplot for the paper
def distchk(p1, p2, xr, yr, w, h, md):
    '''
    Point-to-Point distance calculator
    '''
    sxd = ((float(p2[0]) - float(p1[0])) / xr) * w
    syd = ((float(p2[1]) - float(p1[1])) / yr) * h
    d = math.sqrt((sxd**2+syd**2))
    return d > md

def build_coord_v2(i, xr, yr, w, h, md):
    '''
    Takes array inputs q and i, adjusts them for swarmplotting, and returns xvec and yvec appropriately.
    '''
    xvecset = {r:[] for r in i.keys()}
    yvecset = {r:[] for r in i.keys()}
    #i = y, r = x
    for r in i.keys():
        #index q and i with r to get the arrays that need to be distchk'd against each other.
        #the adjusted set will be added to xvec
        for arrayind in range(0, len(i[r])):
            safe = False
            scale = 0
            while not safe:
                #inital value is checked twice, which is wasteful if its not safe in its initial location but whatever
                #after that the xdir scale values are different.
                for xdir in [-1,1]:
                    check = []
                    xpt = r + scale * xdir
                    for index in range(0, len(xvecset[r])):
                        if abs(i[r][arrayind] - yvecset[r][index]) < md * yr / h: #extra filtering layer
                            check.append(distchk([xpt, i[r][arrayind]], [xvecset[r][index], yvecset[r][index]], xr, yr, w, h, md))
                    if all(check):
                        safe = True 
                    if safe == True: #hypothetically, I shouldn't need this. Yet I do anyways, or it hangs forever.
                        break 
                scale += md * xr / w
            if scale > .5 * xr / w:
                print("Point plotting failed")
                continue
            xvecset[r].append(xpt)
            yvecset[r].append(i[r][arrayind])
    
    return xvecset, yvecset 

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-f', '--file', help = 'Path to the file with common and rare breakpoint information', default = '../common_inv_all.txt')
    parser.add_argument('-x', '--fixed', help = 'Path to the file containing fixed breakpoints (different format from the other file', default = '../fixed_inv_2.txt')       
    args = parser.parse_args()
    return args

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
    dups = tandem_dup(crdata)
    dups = pd.DataFrame(dups)
    dups.columns = ["TanDup"]
    crdata = pd.concat([crdata, dups], 1)
    coms = list(crdata.loc[crdata['Freq'] == 'Common']['TanDup'])
    rares = list(crdata.loc[crdata['Freq'] == 'Rare']['TanDup'])
    fixed = list(crdata.loc[crdata['Freq'] == 'Fixed']['TanDup'].astype(int))

    from scipy.stats import mannwhitneyu
    print('Commons Longer Than Rares:', mannwhitneyu(coms,rares,alternative='greater')[1])
    print('Fixed Longer Than Commons:', mannwhitneyu(coms,fixed,alternative='less')[1])
    print('Fixed Longer Than Rares:', mannwhitneyu(rares,fixed,alternative='less')[1])

    #making a plot for the paper.
    plt.rcParams.update({'font.size': 22})
    plt.figure(figsize = (10,8))
    panel1 = plt.axes([1/6,1/8,4/6,3/4])
    plt.title('Duplication Lengths')
    panel1.set_xlim(.5,3.5)
    # panel1.set_xlabel('Inversion Frequency Type')
    panel1.set_ylabel(r'log$\mathregular{_{10}}$(bp)')
    panel1.set_yticks([10000,20000,30000,40000,50000,60000])
    panel1.set_yticklabels([10,20,30,40,50,60])
    panel1.set_xticks([1,2,3])
    panel1.set_xticklabels(['Rare','Common','Fixed'])
    ycd = {k:[] for k in [1,2,3]}
    for i,ftype in enumerate([rares,coms,fixed]):
        ycd[i+1] = [math.log(f+1,10) for f in ftype]
    sxc, syc = build_coord_v2(ycd, 3, 60000, 1, 2, .001)

    #try logging?
    panel1.set_ylim(0,6)
    panel1.set_yticks([1,2,3,4,5,6])
    panel1.set_yticklabels([1,2,3,4,5,6])

    for r in sxc.keys():
        plt.scatter(x=sxc[r],y=syc[r], facecolor = 'black', linewidth = 0, s = 25, alpha = .5)

    #russ wants boxplots instead of median lines alone
    #boxplots are based on the y-value alone
    #5th percentile, 25th, 50th, 75th, 95th

    for r, values in syc.items():
    #     ymed = st.mean(values)
        order = sorted(values)
        fifth = np.percentile(order, 5)
        twentyfifth = np.percentile(order, 25)
        half = np.percentile(order, 50)
        seventyfifth = np.percentile(order, 75)
        ninetyfifth = np.percentile(order, 95)
        
        panel1.plot([r-1/7,r+1/7],[fifth,fifth], color = 'black', linewidth = 1)
        panel1.plot([r,r],[fifth, twentyfifth], color = 'black', linewidth = 1)
        panel1.plot([r-2/7,r+2/7],[twentyfifth,twentyfifth], color = 'black', linewidth = 1)
        panel1.plot([r-2/7,r-2/7],[twentyfifth,half], color = 'black', linewidth = 1)
        panel1.plot([r+2/7,r+2/7],[twentyfifth,half], color = 'black', linewidth = 1)
        panel1.plot([r-2/7,r+2/7],[half,half], color = 'black', linewidth = 1)
        panel1.plot([r-2/7,r-2/7],[half, seventyfifth], color = 'black', linewidth = 1)
        panel1.plot([r+2/7,r+2/7],[half, seventyfifth], color = 'black', linewidth = 1)
        panel1.plot([r-2/7,r+2/7],[seventyfifth,seventyfifth], color = 'black', linewidth = 1)
        panel1.plot([r,r],[seventyfifth, ninetyfifth], color = 'black', linewidth = 1)
        panel1.plot([r-1/7,r+1/7],[ninetyfifth,ninetyfifth], color = 'black', linewidth = 1)

    plt.savefig('tandup_wlog.png', dpi = 600)

if __name__ == "__main__":
    main()