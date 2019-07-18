#!/usr/bin/python

#because canvas won't take text files, here's the code I used to generate the swarmplots of the figure and the code for the boxplot proper and the figure description below
#################################################################################################################################
# FIGURE DESCRIPTION
#################################################################################################################################

# Figure S1: Insulators associate bidirectionally with inversion breakpoints.
# Panel A shows a swarmplot of distances to the nearest insulator for each breakpoint of that category within a 10kb window. 
# Negative distances are outside the inversion span, positive are inside.
# Grey points are permuted expectation data representing the background; black points are the real values.
# We see no enrichment in any frequency category on either side of the breakpoint, indicating that the insulator association is bidirectional and agnostic to being inside or outside the inversion.
# Panel B shows boxplots for the distribution of insulator distance values with a logarithmic axis.
# We see a strong trend of decreasing insulator distances with increasing population frequency category (rare < common < fixed).

#################################################################################################################################
# CODE FOR SWARMPLOTS
#################################################################################################################################

import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np
import random

def distchk(p1, p2, w, h, md):
    '''
    Point-to-Point distance calculator
    '''
    sxd = (float(p2[0]) - float(p1[0])) / w * h
    syd = (float(p2[1]) - float(p1[1])) / h * w
    d = math.sqrt((sxd**2+syd**2))
    return d > md
        
def scatterset(xvec, md, si, w=9, h=3):
    #takes the set of xvec and piles them up, giving additional ycoord based on minimum distance md
    #base y is 0. For every x, check if its colliding with any other x
    #if it is, add a small bit to the y and check again
#     coords = list(zip([x/1000 for x in xvec], [0 for i in range(0,len(xvec))])) #include conversion to kbp
    coords = [x/1000 for x in xvec]
    xyd = {c:0 for c in coords}
    for i, c in enumerate(coords):
        #start checking if its colliding.
        safe = False
        yshift = 0
        if i > si:
            yshift = .1 #initialize with some height on real (black) points so they're not stuck at the bottom.
        while not safe:
            check = []
            for oi, oc in enumerate(coords):
                if oi != i:
                    check.append(distchk([c, yshift], [oc, xyd[oc]], w, h, md))
#             print(all(check))
            if all(check):
#                 print("Value of {} set: {}".format(i, yshift))
                xyd[c] = yshift
                safe = True
#                 print(len(set(xyd.values())))
                break
            yshift += md / h 
    return [(x,y) for x,y in xyd.items()]
# xys = scatterset(distance_data, .1)

def greypt(classx, bpdata, pnum = 500):
    chromlen = {'2L':23515712, '2R':25288936, '3L':25288936, '3R':32081331, 'X':23544271, '4':1830000} #values as of DM6 assembly, chr 2/3/X only for first testing, chr4 an estimate
    sims = []
    sime = []
    for entry in np.array(bpdata):
        # print('bplen', pnum%len(bpdata))
        # print(int((pnum - pnum%len(bpdata))/len(bpdata)))
        for p in range(0,int(
            (pnum - pnum%len(bpdata)
            )/len(bpdata)
            )):
            nent = [e for e in entry]
            nent[2] = random.choice(range(0,chromlen[entry[0]] - abs(int(entry[3])-int(entry[2]))))
            nent[3] = nent[2] + abs(int(entry[3])-int(entry[2]))
    #             print(entry)
    #             print(nent)
            dset = insulator_track_listv(classx, [nent])
            for sd in dset[0]:
                sims.append(sd) 
            for se in dset[1]:
                sime.append(se)
    return (sims, sime)

def insulator_track_listv(classx, bpdata, cutoff = 20000):
    dvec_s = []
    dvec_e = []
    for spent in np.array(bpdata):
#         print(spent)
        minis = []
        minie = []
        # print(len(classx[spent[0]]))
        for element in classx[spent[0]]:
            # print(spent)
            if abs(int(element) - int(spent[2])) < cutoff: 
                if int(spent[2]) <= int(element) < int(spent[3]):
                    sign = -1
                else:
                    sign = 1
                minis.append((int(element) - int(spent[2]))*sign)
            if abs(int(element) - int(spent[3])) < cutoff:
                if int(spent[2]) <= int(element) < int(spent[3]):
                    sign = -1
                else:
                    sign = 1
                minie.append((int(element)- int(spent[3]))*sign)

#         if abs(minis) < abs(minie)
        try:
            dvec_s.append(minis*sign)
#         else:
            dvec_e.append(minie*sign)
        except:
            continue
            # print(minis, minie)
    
    return dvec_s, dvec_e

def main():
    chromlen = {'2L':23515712, '2R':25288936, '3L':25288936, '3R':32081331, 'X':23544271, '4':1830000} #values as of DM6 assembly, chr 2/3/X only for first testing, chr4 an estimate
    common_rares = []
    with open('common_inv_all.txt', 'r') as cinv:
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
    def tandem_dup(crdata):
        dups = []
        for entry in np.array(crdata):
            dups.append(abs(int(entry[2])-int(entry[3])))
        return dups
    dups = tandem_dup(crdata)
    dups = pd.DataFrame(dups)
    dups.columns = ["TanDup"]
    crdata = pd.concat([crdata, dups], 1, sort = False)
    fixed_singles = []
    with open("fixed_inv_2.txt", 'r') as fixinv:
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
    fixed_singles = fixed_singles.drop(19) #problematic row
    crdata = pd.concat([crdata, fixed_singles], ignore_index = True, sort = False)
    crdata = crdata[['Chro','Label','Forward','Reverse','Freq','TanDup']]
    with open('public_datasets/GSE26905_Dm_Insulator_Classes_r6.GFF3', 'r') as insulators:
        #data structure is two dicts, key of chromosome, value of a set of insulator element points
        #treating as point bases because all binding sites are 10 bases long
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
    for ftype in ['Fixed','Common','Rare']:
        print(ftype)
        bpdata = crdata.loc[crdata['Freq'] == ftype]
        ds, de = insulator_track_listv(class1, bpdata)
        fs, fe = greypt(class1, bpdata)
        print('permuted data generated.')
        d_data = ds
        f_data = fs
        for dset in de:
            dflip = []
            if len(dset)>0:
                dflip = [e*-1 for e in dset]
            d_data.append(dflip)
        for fset in fe:
            fflip = []
            if len(fflip)>0:
                fflip = [e*-1 for e in fset]
            f_data.append(fflip)
        ndat = []
        for dlist in d_data:
            for d in dlist:
                ndat.append(d)
        fdat = []
        for flist in f_data:
            for f in flist:
                fdat.append(f)
        switchind = len(fdat)
        distance_data = fdat + ndat
    #     distance_data = fdat + ndat
        xys = scatterset(distance_data, .1, si = switchind)
        print('offsets calculated')
        plt.figure(figsize = [6,3])
        panel1 = plt.axes([1/12,1/3, .75,.5])
        if ftype == 'Rare':
            panel1.set_title("Short-Range Insulator Distance Distribution")
        panel1.set_ylabel(ftype)
        panel1.set_xlim(-10,10)
        panel1.set_ylim(0,.9)
        panel1.set_yticks([])
        panel1.set_xlabel('Distance [kbp]')
        ##add an extra panel, zoomed into the break
    #     panel2 = plt.axes([10/12,1/3,.2,.5])
    #     panel2.set_xlim(-10,10)
    #     panel2.set_ylim(0,.9)
    #     panel2.set_yticks([])
    #     panel2.set_xlabel('Distance (kbp)')
        panel1.axvline(x=0, color= 'black', linewidth = .1)
    #     panel2.axvline(x=0, color= 'black', linewidth = .1)
        xs = [x[0] for x in xys]
        ys = [y[1] for y in xys]
        print([y for y in ys[switchind:] if y < .05])
        panel1.scatter(xs[:switchind], ys[:switchind], color = 'grey', s = 25)
        panel1.scatter(xs[switchind:], ys[switchind:], color = 'black', s = 25)
        plt.savefig(ftype + "_test.png", dpi = 600)

#################################################################################################################################
# PARTIAL CODE FOR BOXPLOT 
#################################################################################################################################

plt.figure(figsize = (6,10.8))
panel1 = plt.axes([1/6,1/4,2/3,1/2])
plt.title('Distances to Insulators by Frequency')
panel1.set_xlim(.5,3.5)
# panel1.set_xlabel('Inversion Frequency Type')
# panel1.set_ylabel('Distance (log10(kbp))')
panel1.set_ylabel('Distance [bp]')
panel1.set_xticks([1,2,3])
panel1.set_xticklabels(['Rare','Common','Fixed'])
panel1.set_ylim(0,5)
# panel1.set_yticks([10000,20000,30000,40000,50000])
panel1.set_yticks([1,2,3,4,5])
panel1.set_yticklabels(["$10^1$", "$10^2$", "$10^3$", "$10^4$", "$10^5$"])
ycd = {k:[] for k in [1,2,3]}
for i,ftype in enumerate([a,c,x,f]):
    ycd[i+1] = ftype
syc = {k:[] for k in ycd.keys()}
for cat, darray in ycd.items():
    for d in darray:
        syc[cat].append(math.log(d+1, 10))
for r, values in syc.items():
    order = sorted(values)
    fifth = np.percentile(order, 5)
    twentyfifth = np.percentile(order, 25)
    half = np.percentile(order, 50)
    seventyfifth = np.percentile(order, 75)
    ninetyfifth = np.percentile(order, 95)
    print(fifth, twentyfifth, half, seventyfifth, ninetyfifth)
    
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

    
#add a box representing the expected distribution core.    
# stf = st.median([abs(v) for v in f])
# print(stf)
stf = math.log(st.median(f), 10)
panel1.axhline(y=stf, color = (.5,.5,.5), linewidth = 1, linestyle = '--')

plt.savefig('insulator_test.png', dpi = 600)