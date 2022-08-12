#!/usr/bin/env python3
import re
current = {'strand':0,'frame':0,'ovlap':0}
counts = {'strand':[],'frame':[],'ovlap':[]}
prv = {'start':None,'stop':None,'frame':None,'strand':None}
with open(snakemake.input[0], 'r') as f:
    for line in f:
        line = re.sub('<|>', '', line)
        l = line.split()
        coords = l[0].split('..')
        frame = str(int(coords[0]) % 3) + l[1]
        if prv['start'] is None:
            prv['start'] = coords[0]
            prv['stop'] = coords[1]
            prv['strand'] = l[1]
            prv['frame'] = frame
            continue
        if l[1] == prv['strand']:
            current['strand'] += 1
            if frame == prv['frame']:
                current['frame'] += 1
            else:
                counts['frame'].append(current['frame'])
                current['frame'] = 0
        else:
            counts['strand'].append(current['strand'])
            counts['frame'].append(current['frame'])
            current['strand'] = 0
            current['frame'] = 0
        if coords[0] < prv['stop']:
            current['ovlap'] += 1
        else:
            counts['ovlap'].append(current['ovlap'])
        prv['start'] = coords[0]
        prv['stop'] = coords[1]
        prv['strand'] = l[1]
        prv['frame'] = frame
st = open(snakemake.output.strnd,'w')
fr = open(snakemake.output.frame,'w')
ov = open(snakemake.output.ovlps,'w')
st.write('\n'.join([str(x) for x in counts['strand']]))
fr.write('\n'.join([str(x) for x in counts['frame']]))
ov.write('\n'.join([str(x) for x in counts['ovlap']]))
st.close()
fr.close()
ov.close()