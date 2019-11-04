#!/usr/bin/env python3

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
#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-f', '--file', help = 'Path to the file with common and rare breakpoint information', default = 'common_inv_all.txt')
    parser.add_argument('-x', '--fixed', help = 'Path to the file containing fixed breakpoints (different format from the other file', default = 'fixed_inv_2.txt')
    parser.add_argument('-c', '--chromatin', help = 'Path to the chromatin file', default = 'droso_chromatin_r6.GFF3')
    
    args = parser.parse_args()
    return args

#general function and global variables
def tandem_dup(crdata):
    dups = []
    for entry in np.array(crdata):
        dups.append(abs(int(entry[2])-int(entry[3])))
    return dups

chromlen = {'2L':23515712, '2R':25288936, '3L':28110227, '3R':32081331, 'X':23544271, '4':1830000} #values as of DM6 assembly, chr 2/3/X only for first testing, chr4 an estimate

#chromatin functions and classes
def GFFRead(filename):
    '''
    Function that translates GFF3 formatted input files into an array format, with a row for each base.
    '''
    array = []
    with open(filename) as gff:
        for entry in gff.readlines()[5:]:
            ent = entry.strip().split()
            fixentry = [ent[0], ent[3], ent[4], ent[5], ' '.join(ent[8:10])]
            array.append(fixentry)
    return [filename[:-4], array] #return the original name of the file sans the .gff extension


class chromatin:
    '''
    Takes an array of annotated elements in the format chr, start, stop, type, architecture (if entries are strings, start and stop are coerced to int)
    Returns a statistical summary of the elements.
    Used to guide further analysis.
    '''
    def __init__(self, eleset):
        self.chrom = {}
        for entry in eleset:
            try:
                chro = entry[0]
                start = entry[1]
                stop = entry[2]
                ctype = int(entry[3][:-2])
                self.architecture = entry[4]
        #         arch = entry[4] 
                #to make the following code more readable while I'm doing design. Entry 4 not strictly necessary at this point.
                tmp = self.chrom.get(chro, {})
                tmp2 = tmp.get(ctype, [])
                tmp[ctype] = tmp2
                self.chrom[chro] = tmp
                tmp2.append((int(start), int(stop)))
                self.chrom[chro][ctype] = tmp2
            except:
                continue
    
    def lenstat(self, region, ctype, graph = False):
        '''
        Taking a region (2L or whatever) and chromatin type as a string and int input, retrives the relevant set and generates a summary of the elements average length and std deviation of said length.
        '''
        import statistics as st
        cset = self.chrom[region][ctype]
        lengths = []
        for entry in cset:
            lengths.append(entry[1]-entry[0])
        if graph:
            import seaborn as sns
            import matplotlib.pyplot as plt
            plt.figure()
            sns.distplot(lengths)
            plt.show()
        return len(lengths), st.mean(lengths), st.stdev(lengths)


class pixelate:
    '''
    This is the class that will calculate type proportions for a resolution size in a dataset.
    It will include methods to select and graph an area, or return a summary of a whole chromosome.
    '''
    def __init__(self, data, resolution = 10000):
        '''
        Feed in data and a resolution size to generate a pixelate object, which will contain internal objects representing pixel parameters.
        Data is expected in double dictionary format of chromosome-ctype ala chromatin.chrom.
        '''
        self.resolution = resolution
        chromlen = {'2L':23515712, '2R':25288936, '3L':25288936, '3R':32081331, 'X':23544271, '4':1830000} #values as of DM6 assembly, chr 2/3/X only for first testing, chr4 an estimate
        self.pixels = {c:{p:{t:0 for t in range(1,10)} for p in range(0, chromlen[c]+self.resolution, self.resolution)} for c in chromlen.keys()} #a triple nested dictionary - chromosome to pixel start coordinate to chromatin type.
        for chrom, ctype in data.items():
            for cname, ltups in ctype.items():
                for loc in ltups:
                    try:
                        #I want to find the pixel/s that contain the region covered by this loc
                        #calculate how much of the region they own, and add it
                        #using string manipulation, because I'm always rounding down!
                        if loc[1] < loc[0]:
                            loc = [loc[1], loc[0]] #make sure they're smaller to bigger. just to be safe.
                        lpixs = self._conv(loc[0], len(str(resolution))-1)
                        lpixe = self._conv(loc[1], len(str(resolution))-1)
                        if lpixs == lpixe:
                            prev = self.pixels[chrom][lpixs].get(cname, 0)
                            self.pixels[chrom][lpixs][cname] = prev + (loc[1]-loc[0])
                        else:
                            #need to define a range of pixels that are fully included in this chromatin span and add a full complement of 10k bases to each.
                            for pix in range(lpixs+resolution, lpixe):
                                assert pix != lpixs and pix != lpixe
                                prev = self.pixels[chrom][pix].get(cname, 0)
                                self.pixels[chrom][pix][cname] = prev + resolution
                            sprev = self.pixels[chrom][lpixs].get(cname, 0)
                            eprev = self.pixels[chrom][lpixe].get(cname, 0)
                            self.pixels[chrom][lpixs][cname] = sprev + (lpixs + resolution - loc[0]) #the bases before the next pixel, thus in this one.
                            self.pixels[chrom][lpixe][cname] = eprev + (loc[1] - lpixe) #the bases after the start of the second pixel.
                    except KeyError:
#                         print("Location error: {} {} {}".format(chrom, cname, loc))
                        continue
    def _conv(self, loc, sigs = 4):
        '''
        Internal method that uses string manipulation to convert an integer input into a lower entry that would be a valid pixel for resolution.
        '''
        return int(str(loc)[:-sigs] + '0' * sigs)
        
    def graph(self, chrom, start, stop, *ctypes):
        '''
        Takes a start and stop window frame; generates and displays a proportion graph for each called type across the window.
        '''
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns
        pixs = self._conv(start)
        pixe = self._conv(stop)
        props = {k:[] for k in range(1,10)} #key is chromatin type, value is a proportion array.
        for ctype in ctypes:
            for lpix in range(pixs, pixe+self.resolution, self.resolution):
                props[ctype].append(self.pixels[chrom][lpix][ctype])
        data = [[i] + [props[c][i] for c in ctypes] for i in range(len(props[ctypes[0]]))]
        data = pd.DataFrame(data)
        data = pd.melt(data, id_vars = 0, var_name = "ctype", value_name = "prop")
        sns.catplot(x = 0, y = 'prop', hue = 'ctype', data = data, kind = "bar")
        plt.show()
    
    def validate(self):
        '''
        Method calculates pixel heterogeneity and returns the average. 
        A very high value indicates that the resolution set for pixelate may be too high (individual pixels are not capturing entire regions).
        A low value indicates that the resolution set for pixelate may be too low (multiple regions are being captured and homogenized).
        '''
        import statistics as st
        pavg = []
        for chrom, pixset in self.pixels.items():
            for lpix, props in pixset.items():
                #props is a dictionary with keys of chromatin states and values of their counts.
                preset = [v for v in props.values() if v >= 0]
#                 print(props)
                #purity = proportion in the largest vs in the rest.
                if sum(preset) != 0:
                    pavg.append(max(preset)/sum(preset))
        return st.mean(pavg)

    def fetch_pixel(self, chrom, point):
        #get the pixel containing a given point and return its values.
        #this is a dumb/slow algo 
        p2ret = 0
        for pix in self.pixels[chrom].keys():
            try:
                if int(pix) > int(point):
                    return self.pixels[chrom][p2ret]
                else:
                    p2ret = pix
            except:
                continue


def getchrom(bps, pixels):
    '''Iteratively call pixelate.fetch_pixel() on a series of bps locations.'''
    for entry in bps:
        try:
            yield pixels.fetch_pixel(entry[0], (int(entry[2]) + int(entry[3]))/2)
        except:
            try:
                yield pixels.fetch_pixel(entry[0][0], (entry[0][2] + entry[0][3])/2)
            except:
                print(entry)
                
def sumchrom(bps, pixels):
    '''Summarize chromatin across a series of breakpoints based on a pixelate object.'''
    cset = getchrom(bps, pixels)
    vals = {k:[] for k in range(1,10)}
    # print(cset)
    for tdict in cset:
        for ctype, val in tdict.items():
            vals[ctype].append(val)
    for k,v in vals.items():
        try:
            k = [st.mean(v), st.stdev(v)]
        except:
            k = [v, 'N/A']
    return vals

def makedist(bps, pixels):
    '''Calls getchrom and packs the data nicely'''
    import pandas as pd
    cset = getchrom(bps, pixels)
    vals = {k:[] for k in range(1,10)}
    for tdict in cset:
        for ctype in range(1,10):
            tdict[ctype] = tdict.get(ctype, 0)
        for ctype, val in tdict.items():
            vals[ctype].append(val)
    data = pd.DataFrame(vals)
    return data

def chrombinary(ninestate):
    '''
    Take a nine-state model matrix and convert it to a trinary vector of active/inactive/mixed as strings.
    '''
    trinvec = []
    for entry in ninestate:
        #calculate the relative proportion in active states (1-5) and inactive (6-9).
        #>.7 in active = active, >.7 in inactive = inactive, <.7 in both = mixed.
        esum = sum([int(e) for e in entry])
        active = sum([int(e) for e in entry[0:6]])
        inactive = sum([int(e) for e in entry[6:10]])
#         print(active, inactive)/
        assert inactive + active == esum
        if esum == 0:
            trinvec.append("N/A")
        elif active/esum >= .9:
            trinvec.append("Active")
        elif inactive/esum >= .9:
            trinvec.append("Inactive")
        else:
            trinvec.append("Mixed")
    return trinvec

def chromo_permuter_trunc(crdata, pixels, pnum = 5):
    '''
    Another permuter test; compares permuted crdata breakpoints to the set of chromatin in pixels. 
    '''
    perset = pd.DataFrame(columns = crdata.columns)
    for i, entry in crdata.iterrows():
        for p in range(0, pnum):
            nent = entry
            nent['Forward'] = random.choice(range(0,chromlen[entry['Chro']] - abs(int(entry['Forward'])-int(entry['Reverse']))))
            nent['Reverse'] = nent[2] + abs(int(entry['Forward'])-int(entry['Reverse']))
            perset.loc[i*pnum+p] = nent
#     print(perset)
    return makedist(np.array(perset), pixels)

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

    #define a dataframe that will contain the break data for the rest of the analysis
    crdata = pd.DataFrame(common_rares[:-3])
    crdata.columns = ['Chro','Label','Forward','Reverse','Freq']
    gff = GFFRead(args.chromatin)[1]
    nines = []
    thirties = []
    for entry in gff:
        if int(entry[4][3]) == 9:
            nines.append([entry[0][:2]] + entry[1:4] + ['9_state'])
        elif int(entry[4][3]) == 3:
            thirties.append([entry[0][:2]] + entry[1:4] + ['30_state'])
    #analysis proceeds with only nine-state data from this point; uses nine-state to determine active/inactive/mixed status of region.    
    test = chromatin(nines)
    pixels = pixelate(test.chrom, 10000)
    chrodist= makedist(np.array(crdata), pixels)
    crdata = pd.concat([crdata, chrodist], 1)

    trinvec = chrombinary(np.array(crdata.loc[:,range(1,10)]))
    rerep = {k:0 for k in set(trinvec)}
    for v in trinvec:
        rerep[v] += 1
    # print(rerep)
    # print(trinvec)
    trinvec = pd.DataFrame(trinvec)
    trinvec.columns = ["Activity"]
    crdata = pd.concat([crdata, trinvec], 1)
    if args.verbose:
        print("Inversion and Chromatin Data as Read")
        print(crdata)
    arrays = {'freq':[], 'act':[], 'count':[]}
    for ftype in ['Rare','Common','Fixed']:
        bpdat = crdata.loc[crdata['Freq'] == ftype]
        for atype in ['Active','Mixed','Inactive']:
            bpdat2 = bpdat.loc[bpdat['Activity'] == atype]
            arrays['freq'].append(ftype)
            arrays['act'].append(atype)
            arrays['count'].append(len(bpdat2))
    sns.barplot(x='freq',y='count',hue='act',data = arrays)
    plt.show(block=False)

    #check for association of blending (PEV) inversions with frequency
    names = list(set(crdata['Label']))
    pev = {'Fixed':{'PEV':0, 'Non':0}, 'Common':{'PEV':0, 'Non':0}, 'Rare':{'PEV':0, 'Non':0}, 'Fake':{'PEV':0, 'Non':0}}
    for n in names:
        bpdata = crdata.loc[crdata['Label'] == n]
        acts = set(bpdata['Activity'].tolist())
        if len(acts) > 1 and not 'Mixed' in acts:
            pev[bpdata['Freq'].tolist()[0]]['PEV'] += 1
        else:
            pev[bpdata['Freq'].tolist()[0]]['Non'] += 1
        expect_raw = chromo_permuter_trunc(crdata.loc[crdata['Freq'] == bpdata['Freq'].tolist()[0]], pixels, pnum = 10)
        expect = [e for e in chrombinary(np.array(expect_raw)) if e != "N/A"]
        for p in range(0,10):
            fpair = np.random.choice(expect, size = 2, replace = False)
            acts = set(fpair)
            if len(acts) > 1 and not 'Mixed' in acts:
                pev['Fake']['PEV'] += 1
            else:
                pev['Fake']['Non'] += 1

    names = list(set(crdata['Label']))
    cpev = {'Fixed':{'PEV':0, 'Non':0}, 'Common':{'PEV':0, 'Non':0}, 'Rare':{'PEV':0, 'Non':0}}
    pevdist = {'PEV':[], 'Non':[]}
    for n in names:
        bpdata = crdata.loc[crdata['Label'] == n]
        acts = set(bpdata['Activity'].tolist())
        if len(acts) > 1 and not 'Mixed' in acts:
            cpev[bpdata['Freq'].tolist()[0]]['PEV'] += 1
        else:
            cpev[bpdata['Freq'].tolist()[0]]['Non'] += 1
    
    #this part is pretty slow
    for p in range(1000):
    #     print(p)
        pc = {'Fixed':0, 'Common':0, 'Rare':0}
        nc = {'Fixed':0, 'Common':0, 'Rare':0}
        for n in names:
            bpdata = crdata.loc[crdata['Label'] == n]
            expect_raw = chromo_permuter_trunc(crdata.loc[crdata['Freq'] == bpdata['Freq'].tolist()[0]], pixels, pnum = 1)
            expect = [e for e in chrombinary(np.array(expect_raw)) if e != "N/A"]
            fpair = np.random.choice(expect, size = 2, replace = False)
            acts = set(fpair)
            if len(acts) > 1 and not 'Mixed' in acts:
                pc[bpdata['Freq'].tolist()[0]] += 1
            else:
                nc[bpdata['Freq'].tolist()[0]] += 1
        pevdist['PEV'].append(pc)
        pevdist['Non'].append(nc)
    xbox = {k:[] for k in ['Fixed', 'Common', "Rare"]}
    for ftype in xbox:
        for i in range(len(pevdist['PEV'])):
            xbox[ftype].append(pevdist['PEV'][i][ftype]/(pevdist['PEV'][i][ftype] + pevdist['Non'][i][ftype]))
    test = pd.DataFrame(xbox)
    test = pd.melt(test, value_vars = ['Fixed','Common','Rare'])
    real = {"Fixed PEV":cpev['Fixed']['PEV']/sum(cpev['Fixed'].values()), 
#         "Fixed NON":cpev['Fixed']['Non']/sum(cpev['Fixed'].values()), 
        "Common PEV":cpev['Common']['PEV']/sum(cpev['Common'].values()), 
#         'Common Non':cpev['Common']['Non']/sum(cpev['Common'].values()), 
        'Rare PEV':cpev['Rare']['PEV']/sum(cpev['Rare'].values()), 
#         "Rare Non":cpev['Rare']['Non']/sum(cpev['Rare'].values())
       }
    mini = pd.DataFrame({'variable':list(real.keys()), 'value':list(real.values())})
    test.rename(index=str, columns={'variable':'Frequency','value':'Proportion of Inversions Blending'}, inplace = True)
    mini.rename(index=str, columns={'variable':'Frequency','value':'Proportion of Inversions Blending'}, inplace = True)
    sns.boxplot(x = 'Frequency', y = 'Proportion of Inversions Blending', data = test, color = 'white')
    sns.scatterplot(x = 'Frequency', y = 'Proportion of Inversions Blending', data = mini, color = 'red', s = 100)
    plt.show(block=False)

    from scipy.stats import fisher_exact
    for f, pevs in pev.items():
        try:
            cf = cpev[f]
            if f != cf:
                mat = [[pevs['PEV'], pevs['Non']], [cf['PEV'], cf['Non']]]
                print('p-value for frequency category',f,'having nonrandom blending')
                print(fisher_exact(mat)[1])
        except:
            continue

if __name__ == "__main__":
    main()    