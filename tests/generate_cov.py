#!/usr/bin/env python
import argparse
import re
import sys
import os
from lxml import etree


def parse_snakemakes(inp, out, uncov):

    ## Read list of uncovered rules
    uncov_r = set()
    of = open(uncov, 'r')
    for l in of:
        uncov_r.add(l.strip())
    of.close()
    root = etree.Element('coverage')
    root.set('version', '1.0')
    packs = etree.SubElement(root, 'packages')
    pack = etree.SubElement(packs, 'package')
    classes = etree.SubElement(pack, 'classes')

    for smk in inp:
        clas = etree.SubElement(classes, 'class')
        etree.SubElement(clas, 'methods')
        lines_c = etree.SubElement(clas, 'lines')
        clas.set('filename', smk)
        clas.set('name', os.path.basename(smk))
        of = open(smk, 'r')
        ct = 1
        rule_hit = False

        for l in of:
            l = l.strip()
            line = etree.SubElement(lines_c, 'line')
            if l.startswith('rule'):
                rule_name = l.replace('rule ','').replace(':','')
                rule_hit = True
                if rule_name in uncov_r:
                    val = 0
                else:
                    val = 1
                line.set('hits', str(val))
            elif rule_hit:
                line.set('hits', str(val))
            elif not rule_hit:
                line.set('hits', "1")
            line.set('number', str(ct))
            ct = ct + 1
        of.close()
    of = open(out, 'w')
    of.write(etree.tostring(root, pretty_print=True).decode('utf-8'))
    of.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Snakemake cov parser')
    parser.add_argument('--inp', dest='inp', type=str, nargs='+',
                        help='list of smks')
    parser.add_argument('--out', dest='out', type=str,
                        help='output filename.')
    parser.add_argument('--uncov', dest='uncov', type=str,
                        help='uncovered rules.')

    args = parser.parse_args()
    parse_snakemakes(args.inp, args.out, args.uncov)
