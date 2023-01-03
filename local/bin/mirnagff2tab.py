#!/usr/bin/env python
#
# Copyright 2015 Paolo Martini <paolo.cavei@gmail.com>
# Copyright 2015 Enrica Calura <enrica.calura@gmail.com>
# Copyright 2015 Gabriele Sales <gbrsales@gmail.com>

from sys import stdin, exit
from optparse import OptionParser
from os.path import basename
import re

lineno = None
missing_value = None
sep = ' '

#imported from vfork
#from vfork.io.util import safe_rstrip
def safe_rstrip(line):
    return line.rstrip('\r\n')

def format_usage(usage):
    def prefix_length(line):
        length = 0
        while length < len(line) and line[length] in (' ', '\t'):
            length += 1
        return length

    lines = usage.split('\n')
    while len(lines) and len(lines[0].strip()) == 0:
        del lines[0]
    while len(lines) and len(lines[-1].strip()) == 0:
        del lines[-1]

    plen = min(prefix_length(l) for l in lines if len(l.strip()) > 0)
    return '\n'.join(l[plen:] for l in lines)


def ignore_broken_pipe(func):
    try:
        func()
    except IOError as e:
        if e.errno == errno.EPIPE:
            sys.exit(0)
        else:
            raise



def main():
    global missing_value
    assert len(sep) == 1

    options, args = parse_args()
    missing_value = options.missing_value

    if options.list_attributes:
        if len(args) > 0: exit('Unexpected argument number.')
        list_attributes()
    else:
        if len(args) < 1: exit('Unexpected argument number.')
        build_gff(args)

def parse_args():
    parser = OptionParser(usage=format_usage('''
      %prog [OPTIONS] ATTRIBUTE.. <GFF/GTF >TSV

      Transform a GTF/GFF into a tsv with al least one selected ATTRIBUTE as columns.
      The gtf/gff have to be without headers.
    '''))

    parser.add_option('-l', '--list-attributes', action='store_true',
                      default=False, help='list all the available attributes')

    parser.add_option('-m', '--missing-value', default='NA', metavar='STRING',
                      help='write this STRING when an attribute is missing ' \
                           '(default: %default)')

    return parser.parse_args()

def list_attributes():
    attrs = set()
    for toks in iterlines():
        attrs |= tags(toks[8])
    print('\n'.join(sorted(attrs)))

def iterlines():
    global lineno

    for lineno, line in enumerate(stdin, 1):
        tokens = safe_rstrip(line).split('\t')
        if len(tokens) == 9:
            yield tokens
        else:
            die('The input file does not have 9 columns')

def die(reason):
    exit('%s at line %d.' % (reason, lineno))

def tags(str):
    return set(a[0] for a in parse_attrs(str))

def parse_attrs(str, rx=re.compile(r';\s*')):
    parts = rx.split(str.strip())
    parts = (p for p in parts if len(p))
    return (split_attr(p) for p in parts)

def split_attr(str, rx=re.compile('(\S+)\s+(?:"([^"]*)"|(.*))$')):
    m = rx.match(str)
    if m:
        name, fst, snd = m.groups()
        return name, fst or snd
    else:
        pos = str.find(" ")
        if pos == -1: pos = len(str)
        die('Invalid attribute "%s"' % str[:pos])

def build_gff(required_tags):
    output = [None]*len(required_tags)

    for toks in iterlines():
        attrs = dict(parse_attrs(toks[8]))

        for idx, name in enumerate(required_tags):
            output[idx] = attrs.get(name, missing_value)

        print ('\t'.join(toks[:8] + output))


if __name__=='__main__':
    ignore_broken_pipe(main)
