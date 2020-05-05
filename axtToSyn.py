#!/usr/bin/env python3
"""
axtToSyn.py

Generates elongated synteny blocks from a pairwise alignment .axt file.
You can specify a min lastz alignment score to filter .axt file by 
the last filed. The min block length only considers elongated blocks 
larger than this length.

Jake VanCampen
vancampen@ohsu.edu
"""
import argparse
import os

def get_args():
    '''Define and return command line arguments'''
    parser = argparse.ArgumentParser(
        prog='axtToSyn',
        description='Generates synteny blocks from pairwise genome \
        alignment .axt file by block elongation.')
    parser.add_argument('file',
                        help='Relative path to net.axt alignment file.',
                        type=str)
    parser.add_argument('outfile',
                        help='Path to output synteny blocks file \
                              (default: synblocks.csv)',
                        type=str,
                        default="synblocks.csv")
    parser.add_argument('--min-score',
                        type=int,
                        help='Min alignment score of block to be considered for\
                                elongation (defalt: 1e5)',
                        default=2e5,
                        nargs='?'
                        )
    parser.add_argument('--min-blen',
                        type=int,
                        help='Min block len to be considered for\
                                elongation in the first pass (defalt: 1e3)',
                        default=1e3,
                        nargs='?')
    parser.add_argument('--break-len',
                        type=int,
                        help='Set the length of the synteny breakpoint regions that are output (default: 100bp).',
                        default=1e2,
                        nargs='?'),
    parser.add_argument('tname',
                        type=str,
                        help='Target assembly, used to name breakpoints.',
                        )
    parser.add_argument('qname',
                        type=str,
                        help='Query assembly, used to name breakpoints.',
                        )
    return parser.parse_args()


def axtFilter(file, min_score):
    '''
    Generates alignment fields with
    min_score from .axt alignment file.
    '''
    fi=file
    with open(fi, 'r') as f:
        for l in f:
            l=l.strip().split(" ")
            # skip lines without len 9
            if len(l) != 9:
                continue
            # skip score less than min_score 
            elif int(l[8]) < min_score:
                continue
            else:
                yield(l[1:8])

def write_outfile(syn, outname):
    '''
    Write the result of the last pass 
    to a space delimeted output file
    '''
    with open(outname, 'w') as of:
        for b in syn:
            of.write("\t".join(b))
            of.write("\n")

def first_pass(file, min_len, min_score):
    '''
    Make a first pass at elongation of alignment
    blocks from the .axt file. Only keep adjacent blocks
    that have the same target and query genomes and 
    strand. This allows for detection of inversion 
    breakpoints.
    '''
    # block storage for the first pass
    fp=[]
    # start the first block on the first line
    print("Elongating the synblocks on first pass...")

    # iterate the alignment generator
    b=next(axtFilter(file, min_score))

    # init alignment counter
    inb=0
    
    # look for blocks to elongate
    for l in axtFilter(file, min_score):
        # if adjacent blocks have chroms,strand same
        if b[0]==l[0] and b[3]==l[3] and int(b[4])<int(l[4]) and b[6]==l[6]:
            # extend block in chromosome
            b[2]=l[2]
            b[5]=l[5]
            b[6]=l[6]
            inb+=1
        # only add bs greater than min_len
        elif int(b[2])-int(b[1]) > min_len and inb > 2:
            # append a count that represents
            # the number of alignments
            # contributing to a block
            b.append(str(inb))
            fp.append(b)
            b=l
            inb=0
        else:
            # reset block and block counter
            b=l
            inb=0
    return fp


def second_pass(first_pass, min_len):
    '''
    Make a second pass. Elongating blocks generated from 
    the first pass. Using the same adjacency contraints.
    Elongates blocks made of at least three 
    alignments from the first pass, whose length is greater 
    than 500*min_score from the first pass (500kb default).
    '''
    # block storage for the second pass
    sp = []

    # start the first block of the second pass
    s = first_pass[0]
    for fp in first_pass:
        if s[0]==fp[0] and s[3]==fp[3] and int(s[4])<int(fp[4]) and s[6]==fp[6]:
            # elongate block
            s[2]=fp[2]
            s[5]=fp[5]
            # sum alingment counter from current and prev block 
            s[7]=str(int(s[7])+int(fp[7]))

        # at least three contributing blocks
        # len = min_len*500
        elif int(s[7]) > 2 and int(s[2]) - int(s[1]) > min_len*500:
            sp.append(s)
            s=fp
        else:
            # block reset 
            s=fp
    print(sp[1:10])
    return sp


def third_pass(second_pass):
    '''
    Make a thrid pass. Elongating any remaining adjacent blocks 
    from the second pass that have the same target and query chromosome
    and strand.
    '''
    # storage for the thrid pass
    tp=[]
    t=second_pass[0]
    for sp in second_pass:
        if t[0]==sp[0] and t[3]==sp[3] and int(t[4])<int(sp[4]) and t[6]==sp[6]:
            # target end coord
            t[2]=sp[2]
            # query end coord 
            t[5]=sp[5]
            # alignment counter 
            t[7]=str(int(t[7])+int(sp[7]))
            tp.append(t)
        else:
            # append block and reset
            t=sp
    print(tp[1:10])            
    return tp

def write_breakpoints(tp, target, query, bpfile, bklen):
    '''
    Takes the synteny blocks from the third pass
    and write a file that contains the breakpoint 
    regions.
    '''
    # a list of breakpoint regions 
    chr=""
    # chunk list 
    bps=[]
    # breakpoint counter
    nbps=0
    for l in tp:

        # target chromosome number
        tchrom=l[0].replace("chr", "")

        # query chromosome number
        qchrom=l[3].replace("chr", "")

        # increment breakpoint number
        nbps+=1

        # name ends up in fourth field of output breakpoints BED file
        name=f"{target}_{tchrom}_{query}_{qchrom}_"

        # if it's the start of a new chrom
        if l[0] != chr:
            # remove blocks on chromosome ends
            chr=l[0]
            l=[l[0],str(int(int(l[2])-bklen)),l[2],name+f"B{nbps}_L",l[7],l[6]]
            
            # first remove the last if there is one
            # this makes it so we don't count 
            # chromosome ends as breakpoints
            bps = bps[0:-1]

            # then append next first breakpoint
            bps.append(l)
            

        # otherwise in the same chrom
        elif l[0] == chr:

            # construct the right side breakpoint 
            a=[l[0],l[1],str(int(int(l[1])+bklen)),name+f"B{nbps}_R",l[7],l[6]]
            bps.append(a)
            nbps+=1
            
            # construct the left side breakpoint
            b=[l[0],str(int(int(l[2])-bklen)),l[2],name+f"B{nbps}_L",l[7],l[6]]
            bps.append(b)
        else:
            continue

    write_outfile(bps, bpfile)

def main():

    # get command line args
    args=get_args()

    # initiate first pass 
    infile = args.file
    outfile = args.outfile
    min_len = args.min_blen
    min_score = args.min_score

    # assembly versions
    target = args.tname
    query = args.qname

    breakpoint_len = args.break_len

    print(f"Using min alignment score {min_score}...")
    print(f"Using min block length {min_len}")

    # elongate with three pass 
    fp = first_pass(infile, min_len, min_score)
    sp = second_pass(fp, min_len)
    tp = third_pass(sp)

    # write synblocks 
    write_outfile(tp, outfile)

    bkfile=outfile+".breakpoints.bed"
    # write breakpoints
    write_breakpoints(tp, target, query, bkfile, breakpoint_len) 

    # log stats
    nfirst=len(fp)
    nsecond=len(sp)
    nthird=len(tp)

    print(f"{nfirst} synblocks after first pass\n{nsecond} synblocks after second pass\n{nthird} synblocks after third pass.")

if __name__=="__main__":
    main()
