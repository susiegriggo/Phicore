"""
Convert a genbank file to sequences
"""


import argparse


"""
Read a genbank file and do things with it!
"""

import sys
import re
import binascii
import gzip
from Bio import SeqIO
import pandas as pd

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


def is_gzip(gbkf):
    """
    Is the file compressed?
    :param gbkf:
    :return: true if compressed else false
    """

    with open(gbkf, 'rb') as i:
        return binascii.hexlify(i.read(2)) == b'1f8b'


def genbank_seqio(gbkf, verbose=False):
    """
    Get the parser stream
    :param gbkf: genbank file
    :param verbose:
    :return:
    """

    if is_gzip(gbkf):
        handle = gzip.open(gbkf, 'rt')
    else:
        handle = open(gbkf, 'r')
    return SeqIO.parse(handle, "genbank")

def feature_id(seq, feat):
    """
    Choose the appropriate id for the feature
    :param feat: the feature
    :return: the id
    """

    if 'protein_id' in feat.qualifiers:
        return '|'.join(feat.qualifiers['protein_id'])
    elif 'locus_tag' in feat.qualifiers:
        return "|".join(feat.qualifiers['locus_tag'])
    elif 'db_xref' in feat.qualifiers:
        return '|'.join(feat.qualifiers['db_xref'])
    else:
        return seq.id + "." + str(feat.location)

def feat_to_text(feat, qual):
    if qual in feat.qualifiers:
        return " ".join(feat.qualifiers[qual])
    return "-"

def genbank_to_fna(gbkf):
    """
    Parse a genbank file
    :param gbkf: genbank file
    :param verbose:
    :return: a dict of the sequences
    """

    for seq in genbank_seqio(gbkf):
        yield seq.id, seq.seq

def genbank_to_faa(gbkf, complexheader=False, verbose=False):
    """
    Parse a genbank file
    :param gbkf: the genbank file
    :param complexheader: more detail in the header
    :param verbose: more output
    :return: yield the protein id and sequence
    """

    for seq in genbank_seqio(gbkf):
        for feat in seq.features:
            if feat.type != 'CDS':
                continue
            (start, stop, strand) = (feat.location.start.position, feat.location.end.position, feat.strand)
            prtmtd = {
                'EC_number': "",
                'locus_tag': "",
                'note': "",
                'product': "",
                'protein_id': "",
                'ribosomal_slippage': "",
                'transl_table': 11,
                'translation': ""
            }

            cid = feature_id(seq, feat)

            if complexheader:
                loc = f"{start}_{stop}"
                if strand < 0:
                    loc = f"{stop}_{start}"

                cid += f' [{seq.id}] '
                if 'organism' in seq.annotations:
                    cid += f' [{seq.annotations["organism"]}]'
                cid += f' [{seq.id}_{loc}]'
                if 'product' in feat.qualifiers:
                    cid += f' {feat.qualifiers["product"][0]}'
                else:
                    cid += f' [hypothetical protein]'

            if 'translation' in feat.qualifiers:
                yield seq.id, cid, feat.qualifiers['translation'][0]
            else:
                yield seq.id, cid, str(feat.extract(seq).translate().seq)


def genbank_to_functions(gbkf, verbose=False):
    """
    Parse a genbank file
    :param gbkf: the genbank file
    :param verbose: more output
    :return: yield a tple of [protein id, function]
    """
    for seq in genbank_seqio(gbkf):
        for feat in seq.features:
            if feat.type != 'CDS':
                continue

            cid = feature_id(seq, feat)

            prod = "Hypothetical protein"
            if "product" in feat.qualifiers:
                prod = "|".join(feat.qualifiers['product'])

            yield cid, prod


def genbank_to_orfs(gbkf, complexheader=False, verbose=False):
    """
    Parse a genbank file
    :param gbkf:
    :param complexheader:
    :param verbose:
    :return: a dict of the sequences
    """

    for seq in genbank_seqio(gbkf):
        for feat in seq.features:
            if feat.type != 'CDS':
                continue
            (start, stop, strand) = (feat.location.start.position, feat.location.end.position, feat.strand)
            prtmtd = {
                'EC_number': "",
                'locus_tag': "",
                'note': "",
                'product': "",
                'protein_id': "",
                'ribosomal_slippage': "",
                'transl_table': 11,
                'translation': ""
            }

            cid = feature_id(seq, feat)

            if complexheader:
                loc = f"{start}_{stop}"
                if strand < 0:
                    loc = f"{stop}_{start}"

                cid += f' [{seq.id}] '
                if 'organism' in seq.annotations:
                    cid += f' [{seq.annotations["organism"]}]'
                cid += f' [{seq.id}_{loc}]'
                if 'product' in feat.qualifiers:
                    cid += f' {feat.qualifiers["product"][0]}'
                else:
                    cid += f' [hypothetical protein]'

            yield seq.id, cid, str(feat.extract(seq).seq)


def genbank_to_ptt(gbkf, printout=False, verbose=False):
    """
    Convert the genbank file to a table with the same columns of the ptt file
    :param gbkf: the genbank input file
    :param printout: print the table
    :param verbose: more output
    :return: the table
    """

    res = []

    gire = re.compile('GI:(\d+)')
    cogre = re.compile('(COG\S+)')
    for seq in genbank_seqio(gbkf):
        if verbose:
            sys.stderr.write(f"Parsing {seq.id}\n")
        for feat in seq.features:
            if feat.type != "CDS":
                continue

            gi = "-"
            if gire.match(feat_to_text(feat, 'db_xref')):
                gi = gire.match(feat_to_text(feat, 'db_xref'))[1]

            cog = "-"
            if cogre.match(feat_to_text(feat, 'product')):
                cog = cogre.match(feat_to_text(feat, 'product'))[1]

            gene = feat_to_text(feat, 'gene')
            if gene == "-":
                gene = str(feat.location)

            cid = feature_id(seq, feat)

            thisres = [
                f"{feat.location.start}..{feat.location.end}",
                "+" if feat.strand >= 0 else "-",
                (len(feat.location) / 3) - 1,
                gi,
                gene,
                cid,
                cog,
                feat_to_text(feat, 'product')
            ]

            if printout:
                print("\t".join(map(str, thisres)))

            res.append(thisres)

    return res

def genbank_to_phage_finder(gbkf, verbose=False):
    """
    This is a very specific format used by phage_finder (http://phage-finder.sourceforge.net/documentation.htm)
    - contig id, including >
    - length of the contig
    - gene ID
    - start
    - end
    - function

    :param gbkf: the genbank file
    :param verbose: more output
    :return: yields a tple of this data
    """

    for seq in genbank_seqio(gbkf):
        for feat in seq.features:
            if feat.type != 'CDS':
                continue
            cid = feature_id(seq, feat)
            fn = "Hypothetical protein"
            if 'product' in feat.qualifiers:
                fn = feat_to_text(feat, 'product')
            yield [seq.id, len(seq.seq), cid, feat.location.start, feat.location.end, fn]

def genbank_to_pandas(gbkf, mincontiglen, ignorepartials=True, convert_selenocysteine=False, verbose=False):
    """
    This is a bit of a specific format used by phage_boost. its a simple dataframe with a couple of
    additional columns:
        ['contig',
         'id',
         'start',
         'stop',
         'direction',
         'partial',
         'DNAseq',
         'AAseq',
         'header']
    :param ignorepartials: Ignore any gene call with a frameshift (ie. a stop codon in the middle of the sequence)
    :param convert_selenocysteine: PhageBoost crashes with a selenocysteine protein because it is not included in Biopython
    :param gbkf: Genbank file to parse
    :param verbose: more output
    :return: a pandas data frame
    """

    c = 0
    genes = []
    for seq in genbank_seqio(gbkf):
        if len(seq) < mincontiglen:
            sys.stderr.write(f"Skipped {seq.id} because it's length ({len(seq)}) is less than the minimum contig length ({mincontiglen})")
            continue
        for feat in seq.features:
            if feat.type != 'CDS':
                continue

            tid = seq.id + "_" + str(c)
            partial = 0
            # I don't think this is exactly right
            if 'truncated' in feat.qualifiers:
                partial = 1

            dnaseq = str(feat.extract(seq).seq)
            if len(dnaseq) == 0:
                sys.stderr.write(f"The DNA sequence for {feature_id(seq, feat)} was zero, so skipped")
                continue

            # we just do a de novo translation rather than relying on the translation provided
            # in the genbank file that is often wrong
            trans = str(feat.extract(seq).translate().seq)

            while trans.endswith('*'):
                trans = trans[:-1]

            # partial amino acid codes we should ignore
            paa = {'B', 'Z', 'J', 'X', '*'}

            keeporf = True

            if ignorepartials:
                for aa in paa:
                    if aa in trans:
                        sys.stderr.write(f"There is a {aa} in  {feature_id(seq, feat)} so skipped.")
                        keeporf = False

            if not keeporf:
                continue

            if len(trans) == 0:
                sys.stderr.write(f"The translation for {feature_id(seq, feat)} was zero, so skipped")
                continue

            if convert_selenocysteine:
                trans = trans.replace('U', 'C')
            row = [seq.id, c, feat.location.start.position, feat.location.end.position, feat.strand,
                   partial, dnaseq, trans, tid]
            c += 1

            genes.append(row)

    genecalls = pd.DataFrame(genes, columns=['contig', 'id', 'start', 'stop', 'direction', 'partial', 'DNAseq', 'AAseq',
                                             'header'])

    return genecalls



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-g', '--genbank', help='genbank file', required=True)
    parser.add_argument('-c', '--complex', help='complex identifier line', action='store_true')
    parser.add_argument('-a', '--aminoacids', help="output file for the amino acid sequences (.faa will be appended)")
    parser.add_argument('-n', '--nucleotide', help='output file for nucleotide sequence (.fna will be appended)')
    parser.add_argument('-p', '--ptt', help='output file for the ptt protein table')
    parser.add_argument('-o', '--orfs', help='output file for orfs (.orfs will be appended)')
    parser.add_argument('-f', '--functions', help='output file for two column table of [protein id, function]')
    parser.add_argument('--phage_finder', help='make a phage finder file')
    parser.add_argument('--separate', help='separate genbank entries into different files. In this case the ID is prepended to whatever you provide', action='store_true')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    did = False
    if args.nucleotide:
        if args.separate:
            lastid = None
            out = None
            for sid, seq in genbank_to_fna(args.genbank):
                if sid != lastid:
                    if out:
                        out.close()
                    out = open(f"{args.nucleotide}.{sid}.fna", 'w')
                    lastid = sid
                out.write(f">{sid}\n{seq}\n")
            if out:
                out.close()
        else:
            with open(f"{args.nucleotide}.fna", 'w') as out:
                for sid, seq in genbank_to_fna(args.genbank):
                    out.write(f">{sid}\n{seq}\n")
        did = True

    if args.aminoacids:
        if args.separate:
            lastid = None
            out = None
            for seqid, sid, seq in genbank_to_faa(args.genbank, args.complex, args.v):
                if seqid != lastid:
                    if out:
                        out.close()
                    out = open(f"{args.aminoacids}.{seqid}.faa", 'w')
                    lastid = seqid
                out.write(f">{sid}\n{seq}\n")
            if out:
                out.close()
        else:
            with open(f"{args.aminoacids}.faa", 'w') as out:
                for seqid, sid, seq in genbank_to_faa(args.genbank, args.complex, args.v):
                    out.write(f">{sid}\n{seq}\n")
        did = True

    if args.orfs:
        if args.separate:
            lastid = None
            out = None
            for seqid, sid, seq in genbank_to_orfs(args.genbank, args.complex, args.v):
                if seqid != lastid:
                    if out:
                        out.close()
                    out = open(f"{args.orfs}.{seqid}.orfs", 'w')
                    lastid = seqid
                out.write(f">{sid}\n{seq}\n")
            if out:
                out.close()
        else:
            with open(f"{args.orfs}.orfs", 'w') as out:
                for seqid, sid, seq in genbank_to_orfs(args.genbank, args.complex, args.v):
                    out.write(f">{sid}\n{seq}\n")
        did = True


    if args.ptt:
        r = genbank_to_ptt(args.genbank, False, args.v)
        with open(args.ptt, 'w') as out:
            for l in r:
                out.write("\t".join(map(str, l)))
                out.write("\n")
        did = True

    if args.functions:
        with open(args.functions, 'w') as out:
            for pid, prod in genbank_to_functions(args.genbank, args.v):
                out.write(f"{pid}\t{prod}\n")
        did = True

    if args.phage_finder:
        with open(args.phage_finder, 'w') as out:
            for tple in genbank.genbank_to_phage_finder(args.genbank, args.v):
                out.write("\t".join(map(str, tple)) + "\n")
        did = True

    if not did:
        sys.stderr.write("Please provide either a -n, -a, -o, -p, -f output file! (or all)\n")
