import os
import sys
from gsd import hoomd as gsd
import argparse

parser = argparse.ArgumentParser(description='Concatenate gsd files')
parser.add_argument('output_gsd_file', type=str, help='output gsd file')
parser.add_argument('input_gsds', type=str, nargs='+', help='input gsd files')
parser.add_argument('--stride', type=int, default=1, help='stride')
parser.add_argument('--force',action="store_true",default=False,help='overwrite combined file')
args = parser.parse_args()

output_gsd_file = args.output_gsd_file
if os.path.exists(output_gsd_file) and not args.force:
    print("Error: output gsd already existis")
    #print help string
    parser.print_help()
    sys.exit(1)

frame_list = []
for input_gsd in args.input_gsds:
    frame_list.extend(gsd.open(input_gsd,'rb')[1::args.stride])

out_gsd_fh = gsd.open(output_gsd_file,'wb')
out_gsd_fh.extend(frame_list)
out_gsd_fh.close()
