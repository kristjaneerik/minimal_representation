#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

# converts the pos, ref and alt fields of a table to minimal representation

import sys
import argparse
from minimal_representation import get_minimal_representation

parser = argparse.ArgumentParser(description='Produce sample and sample locus files.')
parser.add_argument('-p', '--pos', dest='pos',
    type=int, help='1-based index of column containing POS' )
parser.add_argument('-r', '--ref', dest='ref',
    type=int, help='1-based index of column containing REF' )
parser.add_argument('-a', '--alt', dest='alt',
    type=int, help='1-based index of column containing ALT' )
parser.add_argument('-H', '--has_header', action="store_true")
parser.add_argument('-o', '--outf', nargs='?', type=argparse.FileType('wb'),
                   default=sys.stdout)
parser.add_argument('-i', '--inf', nargs='?', type=argparse.FileType('rb'),
                   default=sys.stdin)
args = parser.parse_args()

# if there is a header line, just print it right back out
if args.has_header:
    header = args.inf.readline()
    args.outf.write(header)

# for other lines, process the POS, REF and ALT fields to minimal representation
for line in args.inf.readlines():
    cols = line.strip().split('\t')
    pos, ref, alt = [cols[i] for i in [args.pos-1, args.ref-1, args.alt-1]]
    newpos, newref, newalt = get_minimal_representation(int(pos),ref,alt)
    cols[args.pos-1] = newpos
    cols[args.ref-1] = newref
    cols[args.alt-1] = newalt
    args.outf.write("\t".join(map(str,cols))+"\n")

if args.inf is not sys.stdin:
    inf.close()

if args.outf is not sys.stdout:
    outf.close()
