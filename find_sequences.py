#!/usr/bin/env python3

#purpose of this script is to iterate through the directories in @trash and fetch out all the chr data and create fasta files of their sequences.

#import
import argparse
import os
import math
import statistics as st
from scipy.stats import percentileofscore

#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    #add args
    parser.add_argument('-b', '--breakpoints', help = 'Path to file containing breakpoints in label/chr/for/rev format.')
    # parser.add_argument('-d', '--directories', nargs = '+', help = 'Directories containing subdirectories chrX/etc further containing .seq files.')
    parser.add_argument('-s', '--spread', type = int, help = 'Set to a value to get sequence that far up and down stream from the breakpoint. Default is 300', default = 5000)
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    #analysis will be repeated, with a separate fasta file generated for each directory under nargs -d. 
    # for direct in args.directories:
        # find_seqs(direct, args.breakpoints, args.spread)
    #cancel above, doing a new type of analysis.
    #start by parsing the new breakpoints file.
    bps = parse_bps(args.breakpoints)
    with open('bp_divergence.txt','w+') as outf:
        #print header to out table
        print('Name\tPercentile\tInversion Divergence\tMost Diverged Strain\tMost Diverged Value\tLeast Diverged Strain\tLeast Diverged Value\tCount of Missing Assemblies', file = outf)
        percentiles = {}
        for bp in bps:
            if args.verbose:
                print("Examining breakpoint {}".format(bp[0]))
            dirl = bp[5]
            midp = math.floor((int(bp[3]) - int(bp[2]))/2 + int(bp[2]))
            #name, chr, loc1, loc2, strain, pop
            mean_divergences = {}
            path = dirl + '/chr' + bp[1]
            indivs = os.listdir(path)
            for seqf in indivs:
                dps = []
                for seqf2 in indivs:
                    if seqf != seqf2:
                        seq1 = get_seq(midp, os.path.join(dirl, 'chr' + bp[1], seqf), args.spread)                    
                        seq2 = get_seq(midp, os.path.join(dirl, 'chr' + bp[1], seqf2), args.spread)
                        dps.append(calculate_divergence(seq1, seq2))
                mean_divergences[seqf.split('_')[0]] = st.mean(dps)
            #pairwise divergences have been calculated for every entry vs every other entry
            #now retreve the strain in bp and calculate its percentile value compared to the list
            #and append to percentile
            rvs = [(k,v) for k,v in mean_divergences.items() if v != 0] #remove 0 entries
            mindiv = min(rvs, key = lambda x:x[1])
            maxdiv = max(rvs, key = lambda x:x[1])
            ambcount = len([v for v in mean_divergences.values() if v == 0]) #count the number that were removed for diagnostics
            rd = mean_divergences[bp[4]]
            if rd == 0:
                print('Breakpoint {} Strain Ambiguous- Skipping'.format(bp[0]))
                print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(bp[0], 'N/A', "AMBIG", maxdiv[0], maxdiv[1], mindiv[0], mindiv[1], ambcount), file = outf)
                # print('Breakpoint {} Strain Ambiguous- Skipping'.format(bp[0]), file = outf)
                continue
            perc = percentileofscore([v for r,v in rvs], rd)
            percentiles[bp[0]] = perc
            # print(bp[0] + "\t" + str(perc), file = outf)
            #print diagnostics as well
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(bp[0], perc, rd, maxdiv[0], maxdiv[1], mindiv[0], mindiv[1], ambcount), file = outf)
            # print('Most Diverged: {} with {}\nLeast Diverged: {} with {}\nNumber Ambiguous {}\nInversion Match Proportion: {}'.format(maxdiv[0], maxdiv[1], mindiv[0], mindiv[1], ambcount, rd), file = outf)
            # print('#######################', file = outf)
            if args.verbose:
                print('Breakpoint {} is in the {} percentile of mean divergences.'.format(bp[0], perc))
                print('Most Diverged: {} with {}\nLeast Diverged: {} with {}\nNumber Ambiguous {}\nInversion Match Proportion: {}'.format(maxdiv[0], maxdiv[1], mindiv[0], mindiv[1], ambcount, rd))

def calculate_divergence(seq1, seq2):
    match = 0
    mismatch = 0
    assert len(seq1) == len(seq2)
    for i, c in enumerate(seq1):
        oc = seq2[i]
        if c != 'N' and oc != 'N':
            if c == oc:
                match += 1
            else:
                mismatch += 1
    if match != 0:
        prop = mismatch / (match + mismatch)
    else:
        return 0
    return prop


def find_seqs(dirl, bps, spread = 300):
    '''
    Top level function to seek() and extract sequences.
    '''
    nbps = parse_bps(bps) #label is 0, chro is 1, forw is 2, backw is 3
    for bp in nbps:
        #get the path to the directory containing the right arm for any given bp
        path = dirl + '/chr' + bp[1]
        #base the location on the midpoint between the forward and reverse breaks for arbitrary reasons
        midp = math.floor((int(bp[3]) - int(bp[2]))/2 + int(bp[2]))
        #generate a name for the custom fasta where all of these outputs will go, based on the bp label + the directory name
        with open(dirl + '.fa','a+') as outf:
            for seqf in os.listdir(path):
                seqstr = get_seq(midp, os.path.join(dirl, 'chr' + bp[1], seqf), spread)
                label = seqf.split('_')[0] + '_' + bp[0] + '_' + dirl
                write_fasta(seqstr, label, outf)

def get_seq(loc, infile, spread = 300):
    '''
    Seek to the target point in the file and collect sequence to either side based on spread offset.
    '''
    with open(infile, 'r') as inf:
        #seek magic. seq file is simple, no characters to alter count.
        inf.seek(loc-spread)
        seq = inf.read(spread*2)
        return seq

def write_fasta(seqstr, label, outf):
    '''
    Add the sequence entry to the fasta in outf. Outf should be an already-opened object.
    '''
    print('>' + label, file = outf)
    print(seqstr, file = outf)

def parse_bps(bps):
    '''
    Parse breakpoint file back to tuples.
    '''
    bpsl = []
    for entry in open(bps, 'r'):
        spent = entry.strip().split()
        bpsl.append(tuple(spent))
    return bpsl

if __name__ == "__main__":
    main()