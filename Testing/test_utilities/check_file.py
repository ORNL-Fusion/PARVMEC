#!/usr/bin/python

from optparse import OptionParser
import os.path

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="File name to check.", metavar="FILE")
(option, args) = parser.parse_args()

if not os.path.exists(option.filename):
    exit(1)
