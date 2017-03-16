#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from os.path import basename, isfile, abspath, join
import sys

def test():
   pass

if __name__ == "__main__":

   SCRIPT = basename(__file__)

   infile, value, outfile = '', 0,''
   try:
       nucleotide_alignment = sys.argv[1]
       DB_to_prefix = sys.argv[2]
   except:
       sys.exit("USAGE: ./%s %s %s" % (SCRIPT, "cds_alignment", "prefix_to_DB_file"))

   test()
