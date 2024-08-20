#!/usr/bin/env python3
import numpy as np
import argparse
from psihdf import rdh5, wrh5

# Get rid of warning when nan or zero value in divisor
np.seterr(divide='ignore', invalid='ignore')

def argParsing():

    parser = argparse.ArgumentParser(description='hipft_map_diff: Compute comparison metrics between two HipFT H5 map files (including any/all realizations).  The metrics are:   RMSD   MAXABS   MCV(RMSD)  MAPE  MAXAPE')

    parser.add_argument('file1',
                        type=str,
                        help="First h5 file to compare.")

    parser.add_argument('file2',
                        type=str,
                        help="Second h5 to compare (reference).")

    parser.add_argument('-diff',
                        type=str,
                        help="Create diff h5 file of the two input files.")

    parser.add_argument('-dtn',
                        type=str,
                        help="Values of theta above this will be ignored. Can be float or division i.e. [np.pi/100]")

    parser.add_argument('-dts',
                        type=str,
                        help="Values of theta below this will be ignored. Can be float or division i.e. [np.pi/100]")

    parser.add_argument('-idtn',
                        type=int,
                        help="Values of theta above this will be ignored.")

    parser.add_argument('-idts',
                        type=int,
                        help="Values of theta below this will be ignored.")

    return parser.parse_args()


# Beginning of comparison functions ---------------------------------------------------------

def getdiff_rmsd(x,y):
    N = np.size(x)
    d = np.sqrt(np.sum((x - y) ** 2)/N)
    return d

def getdiff_maxabs(x, y):
    d = np.max(np.abs(x - y))
    return d

def getdiff_mcvrmsd(x, y):
    top = np.sum((x - y) ** 2)
    bot = np.sum(y ** 2)
    if bot == 0:
        d = -1
    else:
        d = np.sqrt(top/bot)
    return d

def getdiff_mape(x, y):
    N = np.size(x)
    d = np.sum(np.abs((x - y)/y))/N
    if d == float('inf') or np.isnan(d):
        d = -1
    return d

def getdiff_maxape(x, y):
    N = np.size(x)
    d = np.max(np.abs((x - y)/y))
    if d == float('inf') or np.isnan(d):
        d = -1
    return d


# End of comparison functions ---------------------------------------------------------------

args = argParsing()

t1,t2,t3,fieldData1 = rdh5(args.file1)
_,_,_,fieldData2    = rdh5(args.file2)

fieldData1=np.array(fieldData1)
fieldData2=np.array(fieldData2)

ibot = 0;
itop = -1;

if (args.idts):
    ibot = args.idts
    
if (args.idtn):
    itop = len(t2)-args.idtn

if (args.dts):
    dts = eval(args.dts)
    ibot = (np.abs(t2 - dts)).argmin()
    if t2[ibot] < dts:
        ibot+=1
    
if (args.dtn):
    dtn = eval(args.dtn)
    itop = (np.abs(t2 - dtn)).argmin()
    if t2[itop] > dtn:
        itop-=1

if (ibot != 0 or itop != -1):
    fieldData1 = fieldData1[ibot:itop+1,:]
    fieldData2 = fieldData2[ibot:itop+1,:]
    t2 = t2[ibot:itop+1]

listOfdiff = []
listOfdiff.append(getdiff_rmsd(fieldData1, fieldData2))
listOfdiff.append(getdiff_maxabs(fieldData1, fieldData2))
listOfdiff.append(getdiff_mcvrmsd(fieldData1, fieldData2))
listOfdiff.append(getdiff_mape(fieldData1, fieldData2))
listOfdiff.append(getdiff_maxape(fieldData1, fieldData2))

if (args.diff is not None):
    fieldDataDiff = fieldData1 - fieldData2
    if not args.diff.endswith('.h5'):
        diff_file = args.diff+'.h5'
    else:
        diff_file = args.diff

    wrh5(diff_file,t1,t2,t3,fieldDataDiff)

out = ""

for value in listOfdiff:
    if value != -1 and  value != 0:
        out += "%23.17e," % value
    else: 
        out += "%23d," % value

print(out)
