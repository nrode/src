#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: test demultiplex.py
# Copyright (C) 2014 Institut National de la Recherche Agronomique
# License: GPL-3+
# Author: Timothée Flutre
# Versioning: https://github.com/timflutre/quantgen

# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import sys
import os
import getopt
import time
import datetime
from subprocess import Popen, PIPE, check_output
import math
import gzip
import tempfile
import shutil
import itertools

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Data.IUPACData import ambiguous_dna_values

if sys.version_info[0] == 2:
    if sys.version_info[1] < 7:
        msg = "ERROR: Python should be in version 2.7 or higher"
        sys.stderr.write("%s\n\n" % msg)
        sys.exit(1)
        
progVersion = "1.3.0" # http://semver.org/


class TestDemultiplex(object): #Define a class
    
    def __init__(self): #Define attributes
        self.verbose = 0
        self.pathToProg = ""
        self.testsToRun = ["1", "2", "4a", "4b", "4c","4d"]
        self.clean = True
        
        
    def help(self): #Define a method for the help
        """
        Display the help on stdout.
        
        The format complies with help2man (http://www.gnu.org/s/help2man)
        """
        msg = "`%s' tests demultiplex.py.\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Options:\n"
        msg += "  -h, --help\tdisplay the help and exit\n"
        msg += "  -V, --version\toutput version information and exit\n"
        msg += "  -v, --verbose\tverbosity level (default=0/1/2/3)\n"
        msg += "  -p, --p2p\tfull path to the program to be tested\n"
        msg += "  -t, --test\tidentifiers of test(s) to run (default=1-2-4a-4b-4c-4d)\n"
        msg += "  -n, --noclean\tkeep temporary directory with all files\n"
        msg += "\n"
        msg += "Examples:\n"
        msg += "  %s -p ~/src/demultiplex.py -n\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Report bugs to <timothee.flutre@supagro.inra.fr>."
        print(msg); sys.stdout.flush()
        
        
    def version(self): #Define a method for the version
        """
        Display version and license information on stdout.
        """
        msg = "%s %s\n" % (os.path.basename(sys.argv[0]), progVersion)
        msg += "\n"
        msg += "Copyright (C) 2014 Institut National de la Recherche Agronomique (INRA).\n"
        msg += "License GPL-3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        msg += "\n"
        msg += "Written by Timothée Flutre."
        print(msg.encode("utf8")); sys.stdout.flush()
        
        
    def setAttributesFromCmdLine(self): #Define a method to read 
        """
        Parse the command-line arguments.
        """
        try:
            opts, args = getopt.getopt( sys.argv[1:], "hVv:p:t:n",
                                        ["help", "version", "verbose=",
                                         "p2p=", "test=", "noclean"])
        except getopt.GetoptError as err:
            sys.stderr.write("%s\n\n" % str(err))
            self.help()
            sys.exit(2)
        for o, a in opts:
            if o == "-h" or o == "--help":
                self.help()
                sys.exit(0)
            elif o == "-V" or o == "--version":
                self.version()
                sys.exit(0)
            elif o == "-v" or o == "--verbose":
                self.verbose = int(a)
            elif o == "-p" or o == "--p2p":
                 self.pathToProg = a
            elif o == "-t" or o == "--test":
                self.testsToRun = a.split("-")
            elif o == "-n" or o == "--noclean":
                self.clean = False
            else:
                assert False, "invalid option"
                
                
    def checkAttributes(self): #Check the values given as arguments
        """
        Check the values of the command-line parameters.
        """
        if self.pathToProg == "":
            msg = "ERROR: missing compulsory option --p2p"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not os.path.exists(self.pathToProg):
            msg = "ERROR: can't find file %s" % self.pathToProg
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        for t in self.testsToRun:
            if t not in ["1", "2", "4a", "4b", "4c", "4d"]:
                msg = "ERROR: unknown --test %s" % t
                sys.stderr.write("%s\n\n" % msg)
                self.help()
                sys.exit(1)
                
                
    #==========================================================================
            
            
    def beforeTest(self):
        """
        Create a temporary folder with a unique ID
        """
        cwd = os.getcwd()
        testDir = tempfile.mkdtemp(dir=cwd, prefix="tmp_test_")
        os.chdir(testDir)
        if self.verbose > 0:
            print("temp dir: %s" % os.getcwd()); sys.stdout.flush()
        return cwd, testDir
        
        
    def writeTagFile(self, it):
        with open(it, "w") as itHandle:
            txt = ">ind1\nAAA\n"
            txt += ">ind2\nTTT\n"
            txt += ">ind3\nGGG\n"
            txt += ">ind4\nCCC\n"
            itHandle.write(txt)
            
            
    def launchProg(self, ifq1, ifq2, it, met, dist, re, nci):
        """
        Launch demultiplex.py with tha 'args' arguments
        """
        args = [self.pathToProg, #Path to demultiplex.py program
                "--idir", "./",
                "--ifq1", ifq1,
                "--ifq2", ifq2,
                "--it", it,
                "--ofqp", "test",
                "--met", met,
                "--dist", str(dist),
                "-v", str(self.verbose - 1)]
        if re != "":
            args.append("--re") #Add the restriction enzyme name to the 'args' list only if provided
            args.append(re)
        if nci: #Check if nci is true (=do not cut the tag) and add "--nci" to the 'args' list
            args.append("--nci")
        if self.verbose > 0:
            print(" ".join(args))
        # msgs = check_output(args)
        msgs = Popen(args, stdout=PIPE, stderr=PIPE).communicate()
        if self.verbose > 1:
            if msgs[0] != "":
                print("stdout:")
                print(msgs[0][:-1])
            if msgs[1] != "":
                print("sterr:")
                print(msgs[1][:-1])
        return msgs
        
        
    def afterTest(self, cwd, testDir): #Remove testDir and everything in it
        os.chdir(cwd)
        if self.clean:
            shutil.rmtree(testDir)
            
            
    #==========================================================================
    
    
    def test_met1_prepare(self):
        ifq1 = "reads_R1.fastq.gz"
        ifq2 = "reads_R2.fastq.gz"
        ifq1Handle = gzip.open(ifq1, "w")
        ifq2Handle = gzip.open(ifq2, "w")
        
        # pair 1: both reads have perfect tag of ind 2 at bp 1
        txt = "@INST1:1:FLOW1:2:2104:15343:197391 1:N:0\n"
        txt += "TTT" # tag
        txt += "TCAACCTGGAGTTCCAC\n" # insert
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~\n"
        ifq1Handle.write(txt)
        txt = "@INST1:1:FLOW1:2:2104:15343:197391 2:N:0\n"
        txt += "TTT" # tag
        txt += "GTAGCTGAGATCGGAAG\n" # insert
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~\n"
        ifq2Handle.write(txt)
        
        # pair 2: only read 1 has perfect tag of ind 1 at bp 1
        txt = "@INST1:1:FLOW1:2:2104:15343:197392 1:N:0\n"
        txt += "AAA" # tag
        txt += "TCAACCTGGAGTTCCAC\n" # insert
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~\n"
        ifq1Handle.write(txt)
        txt = "@INST1:1:FLOW1:2:2104:15343:197392 2:N:0\n"
        txt += "ATA" # tag
        txt += "GTAGCTGAGATCGGAAG\n" # insert
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~\n"
        ifq2Handle.write(txt)
        
        ifq1Handle.close()
        ifq2Handle.close()
        
        for f in ["test_ind2_R1.fastq.gz", "test_ind2_R2.fastq.gz",
                  "test_unassigned_R1.fastq.gz", "test_unassigned_R2.fastq.gz"]:
            if os.path.isfile(f):
                os.remove(f)
                
        it = "tags.fa"
        self.writeTagFile(it)
        
        return ifq1, ifq2, it
        
        
    def test_met1_comp(self, msgs):
        if not os.path.exists("test_ind2_R1.fastq.gz") or not os.path.exists("test_ind2_R2.fastq.gz"):
            print("test_met1: fail (1)")
            return
        else:
            with gzip.open("test_ind2_R1.fastq.gz") as inFqHandle1, \
                 gzip.open("test_ind2_R2.fastq.gz") as inFqHandle2:
                reads1 = SeqIO.parse(inFqHandle1, "fastq",
                                     alphabet=IUPAC.ambiguous_dna)
                reads2 = SeqIO.parse(inFqHandle2, "fastq",
                                     alphabet=IUPAC.ambiguous_dna)
                for (read1, read2) in itertools.izip(reads1, reads2):
                    if read1.id != "INST1:1:FLOW1:2:2104:15343:197391" and \
                       read2.id != "INST1:1:FLOW1:2:2104:15343:197391":
                        print("test_met1: fail (2)")
                        return
                    if str(read1.seq) != "TCAACCTGGAGTTCCAC" or \
                       str(read2.seq) != "GTAGCTGAGATCGGAAG":
                        print("test_met1: fail (3)")
                        return
        if os.path.exists("test_ind1_R1.fastq.gz") or \
           os.path.exists("test_ind1_R2.fastq.gz"):
            print("test_met1: fail (3)")
            return
        print("test_met1: pass")
        
        
    def test_met1(self):
        if self.verbose > 0:
            print("launch test met1 ...")
            sys.stdout.flush()
        cwd, testDir = self.beforeTest()
        ifq1, ifq2, it = self.test_met1_prepare()
        msgs = self.launchProg(ifq1, ifq2, it, "1", 0, "", False)
        self.test_met1_comp(msgs)
        self.afterTest(cwd, testDir)
        
        
    #==========================================================================

    

    def test_met2_comp(self, msgs):            # checks that the 2 files for each of the 2 indiv exist
        if not os.path.exists("test_ind2_R1.fastq.gz") or \
           not os.path.exists("test_ind2_R2.fastq.gz") or \
           not os.path.exists("test_ind1_R1.fastq.gz") or \
           not os.path.exists("test_ind1_R2.fastq.gz"):
            print("test_met2: fail (1)")
            return
        else:
            with gzip.open("test_ind2_R1.fastq.gz") as inFqHandle1, gzip.open("test_ind2_R2.fastq.gz") as inFqHandle2:     # Open a file from gzip for indiv 2 (tag TTT)
                l1 = list(SeqIO.parse(inFqHandle1, "fastq",             # For Forward, creates a list with Ind ID and sequence in Fasta or Fastaq
                                      alphabet=IUPAC.ambiguous_dna))
                l2 = list(SeqIO.parse(inFqHandle2, "fastq",             # For Reverse, creates a list with Ind ID and sequence in Fasta or Fastaq
                                      alphabet=IUPAC.ambiguous_dna))
                if len(l1) != 1 or len(l2) != 1:                        #Only one element in each list
                    print("test_met2: fail (2)")
                    return
                if l1[0].id != "INST1:1:FLOW1:2:2104:15343:197391" or \
                   l2[0].id != "INST1:1:FLOW1:2:2104:15343:197391":     # check that sequence names are the same
                    print("test_met2: fail (3)")
                    return
                if str(l1[0].seq) != "TCAACCTGGAGTTCCAC" or \
                   str(l2[0].seq) != "GTAGCTGAGATCGGAAG":           # check that sequences cut as expected
                    print("test_met1: fail (4)")
                    return
            with gzip.open("test_ind1_R1.fastq.gz") as inFqHandle1, \
                 gzip.open("test_ind1_R2.fastq.gz") as inFqHandle2:      # for indiv 1 (tag AAA)
                l1 = list(SeqIO.parse(inFqHandle1, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
                l2 = list(SeqIO.parse(inFqHandle2, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
                if len(l1) != 1 or len(l2) != 1:
                    print("test_met2: fail (5)")
                    return
                if l1[0].id != "INST1:1:FLOW1:2:2104:15343:197392" or \
                   l2[0].id != "INST1:1:FLOW1:2:2104:15343:197392":
                    print("test_met2: fail (6)")
                    return
                if str(l1[0].seq) != "TCAACCTGGAGTTCCAC" or \
                   str(l2[0].seq) != "ATAGTAGCTGAGATCGGAAG":
                    print("test_met1: fail (7)")
                    return
        print("test_met2: pass")
        
        
    def test_met2(self):
        if self.verbose > 0:
            print("launch test met2 ...")
            sys.stdout.flush()
        cwd, testDir = self.beforeTest()
        ifq1, ifq2, it = self.test_met1_prepare()
        msgs = self.launchProg(ifq1, ifq2, it, "2", 0, "", False)
        self.test_met2_comp(msgs)
        self.afterTest(cwd, testDir)
        
        
    #==========================================================================
    
    
    def test_met4_prepare(self):
        ifq1 = "reads_R1.fastq.gz"
        ifq2 = "reads_R2.fastq.gz"
        ifq1Handle = gzip.open(ifq1, "w")
        ifq2Handle = gzip.open(ifq2, "w")
        
        # pair 1: read 1 has perfect tag of ind 2 at bp 1; read 2 has no tag
        txt = "@INST1:1:FLOW1:2:2104:15343:197391 1:N:0\n"
        txt += "TTT" # tag
        txt += "TCAACCTGGAGTTCCAC\n" # insert
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~\n"
        ifq1Handle.write(txt)
        txt = "@INST1:1:FLOW1:2:2104:15343:197391 2:N:0\n"
        txt += ""
        txt += "GTAGCTGAGATCGGAAG\n" # insert
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~\n"
        ifq2Handle.write(txt)
        
        # pair 2: read 1 has perfect tag of ind 1 but at bp 2; read 2 has no tag
        txt = "@INST1:1:FLOW1:2:2104:15343:197392 1:N:0\n"
        txt += "TAAA" # tag with a 1-bp shift
        txt += "TAGCTACATHACTACAT\n" # insert
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~~\n"
        ifq1Handle.write(txt)
        txt = "@INST1:1:FLOW1:2:2104:15343:197392 2:N:0\n"
        txt += ""
        txt += "CTCAGCTGGACTCGACT\n" # insert
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~\n"
        ifq2Handle.write(txt)
        
        ifq1Handle.close()
        ifq2Handle.close()
        
        for f in ["test_ind2_R1.fastq.gz", "test_ind2_R2.fastq.gz",
                  "test_unassigned_R1.fastq.gz", "test_unassigned_R2.fastq.gz"]:
            if os.path.isfile(f):
                os.remove(f)
                
        it = "tags.fa"
        self.writeTagFile(it)
        
        return ifq1, ifq2, it
        
        
    def test_met4a_comp(self, msgs):
        if not os.path.exists("test_ind2_R1.fastq.gz") or \
           not os.path.exists("test_ind2_R2.fastq.gz"):
            print("test_met4a: fail (1)")
            return
        else:
            with gzip.open("test_ind2_R1.fastq.gz") as inFqHandle1, \
                 gzip.open("test_ind2_R2.fastq.gz") as inFqHandle2:
                l1 = list(SeqIO.parse(inFqHandle1, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
                l2 = list(SeqIO.parse(inFqHandle2, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
                if len(l1) != 1 or len(l2) != 1:
                    print("test_met4a: fail (2)")
                    return
                if l1[0].id != "INST1:1:FLOW1:2:2104:15343:197391" or \
                   l2[0].id != "INST1:1:FLOW1:2:2104:15343:197391":
                    print("test_met4a: fail (3)")
                    return
                if str(l1[0].seq) != "TCAACCTGGAGTTCCAC" or \
                   str(l2[0].seq) != "GTAGCTGAGATCGGAAG":
                    print("test_met4a: fail (4)")
                    return
        if os.path.exists("test_ind1_R1.fastq.gz") or \
           os.path.exists("test_ind1_R2.fastq.gz"):
            print("test_met4a: fail (5)")
            return
        print("test_met4a: pass")
        
        
    def test_met4a(self):
        if self.verbose > 0:
            print("launch test met4a ...")
            sys.stdout.flush()
        cwd, testDir = self.beforeTest()
        ifq1, ifq2, it = self.test_met4_prepare()
        msgs = self.launchProg(ifq1, ifq2, it, "4a", 0, "", False)
        self.test_met4a_comp(msgs)
        self.afterTest(cwd, testDir)
        
        
    #==========================================================================
    
    
    def test_met4b_comp(self, msgs):
        if not os.path.exists("test_ind2_R1.fastq.gz") or \
           not os.path.exists("test_ind2_R2.fastq.gz"):
            print("test_met4b: fail (1)")
            return
        else:
            with gzip.open("test_ind2_R1.fastq.gz") as inFqHandle1, \
                 gzip.open("test_ind2_R2.fastq.gz") as inFqHandle2:
                l1 = list(SeqIO.parse(inFqHandle1, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
                l2 = list(SeqIO.parse(inFqHandle2, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
                if len(l1) != 1 or len(l2) != 1:
                    print("test_met4b: fail (2)")
                    return
                if l1[0].id != "INST1:1:FLOW1:2:2104:15343:197391" or \
                   l2[0].id != "INST1:1:FLOW1:2:2104:15343:197391":
                    print("test_met4b: fail (3)")
                    return
                if str(l1[0].seq) != "TCAACCTGGAGTTCCAC" or \
                   str(l2[0].seq) != "GTAGCTGAGATCGGAAG":
                    print("test_met4b: fail (4)")
                    return
        if not os.path.exists("test_ind1_R1.fastq.gz") or \
           not os.path.exists("test_ind1_R2.fastq.gz"):
            print("test_met4b: fail (5)")
            return
        else:
            with gzip.open("test_ind1_R1.fastq.gz") as inFqHandle1, \
                 gzip.open("test_ind1_R2.fastq.gz") as inFqHandle2:
                l1 = list(SeqIO.parse(inFqHandle1, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
                l2 = list(SeqIO.parse(inFqHandle2, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
                if len(l1) != 1 or len(l2) != 1:
                    print("test_met4b: fail (6)")
                    return
                if l1[0].id != "INST1:1:FLOW1:2:2104:15343:197392" or \
                   l2[0].id != "INST1:1:FLOW1:2:2104:15343:197392":
                    print("test_met4b: fail (7)")
                    return
                if str(l1[0].seq) != "TAGCTACATHACTACAT" or \
                   str(l2[0].seq) != "CTCAGCTGGACTCGACT":
                    print("test_met4b: fail (8)")
                    return
        print("test_met4b: pass")
        
        
    def test_met4b(self):
        if self.verbose > 0:
            print("launch test met4b ...")
            sys.stdout.flush()
        cwd, testDir = self.beforeTest()
        ifq1, ifq2, it = self.test_met4_prepare()
        msgs = self.launchProg(ifq1, ifq2, it, "4b", 10, "", False)
        self.test_met4b_comp(msgs)
        self.afterTest(cwd, testDir)
        
        
    #==========================================================================
    
    
    def test_met4c_prepare(self):
        ifq1 = "reads_R1.fastq.gz"
        ifq2 = "reads_R2.fastq.gz"
        ifq1Handle = gzip.open(ifq1, "w")
        ifq2Handle = gzip.open(ifq2, "w")
        
        # pair 1: read 1 has perfect tag and cut site of ind 2 at bp 2; read 2 has no tag
        txt = "@INST1:1:FLOW1:2:2104:15343:197391 1:N:0\n"
        txt += "ATTT" # tag with a 1-bp shift
        txt += "CAGCCCTGGAGTTCCAC\n" # insert with ApekI cut site
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~~\n"
        ifq1Handle.write(txt)
        txt = "@INST1:1:FLOW1:2:2104:15343:197391 2:N:0\n"
        txt += ""
        txt += "GTAGCTGAGATCGGAAG\n" # insert
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~\n"
        ifq2Handle.write(txt)
        
        # pair 2: read 1 has perfect tag of ind 1 at bp 2 but no cut site; read 2 has no tag
        txt = "@INST1:1:FLOW1:2:2104:15343:197392 1:N:0\n"
        txt += "TAAA" # tag with a 1-bp shift
        txt += "TAGCTACATHACTACAT\n" # insert with no ApekI cut site
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~~\n"
        ifq1Handle.write(txt)
        txt = "@INST1:1:FLOW1:2:2104:15343:197392 2:N:0\n"
        txt += ""
        txt += "CTCAGCTGGACTCGACT\n" # insert
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~\n"
        ifq2Handle.write(txt)
        
        ifq1Handle.close()
        ifq2Handle.close()
        
        for f in ["test_ind2_R1.fastq.gz", "test_ind2_R2.fastq.gz",
                  "test_unassigned_R1.fastq.gz", "test_unassigned_R2.fastq.gz"]:
            if os.path.isfile(f):
                os.remove(f)
                
        it = "tags.fa"
        self.writeTagFile(it)
        
        return ifq1, ifq2, it
        
        
    def test_met4c_comp(self, msgs):
        if not os.path.exists("test_ind2_R1.fastq.gz") or \
           not os.path.exists("test_ind2_R2.fastq.gz"):
            print("test_met4c: fail (1)")
            return
        else:
            with gzip.open("test_ind2_R1.fastq.gz") as inFqHandle1, \
                 gzip.open("test_ind2_R2.fastq.gz") as inFqHandle2:
                l1 = list(SeqIO.parse(inFqHandle1, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
                l2 = list(SeqIO.parse(inFqHandle2, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
                if len(l1) != 1 or len(l2) != 1:
                    print("test_met4c: fail (2)")
                    return
                if l1[0].id != "INST1:1:FLOW1:2:2104:15343:197391" or \
                   l2[0].id != "INST1:1:FLOW1:2:2104:15343:197391":
                    print("test_met4c: fail (3)")
                    return
                if str(l1[0].seq) != "CAGCCCTGGAGTTCCAC" or \
                   str(l2[0].seq) != "GTAGCTGAGATCGGAAG":
                    print("test_met4c: fail (4)")
                    return
        if os.path.exists("test_ind1_R1.fastq.gz") or \
           os.path.exists("test_ind1_R2.fastq.gz"):
            print("test_met4c: fail (5)")
            return
        print("test_met4c: pass")
        
        
    def test_met4c(self):
        if self.verbose > 0:
            print("launch test met4c ...")
            sys.stdout.flush()
        cwd, testDir = self.beforeTest()
        ifq1, ifq2, it = self.test_met4c_prepare()
        msgs = self.launchProg(ifq1, ifq2, it, "4c", 10, "ApeKI", False)
        self.test_met4c_comp(msgs)
        self.afterTest(cwd, testDir)
        
        
  
     #==========================================================================
     #New method to assign a pair as true only if the pattern (tag+cut site) is found in the forward and/or the reverse read
     #Allows a shift (='dist') in the start of the pattern       
    
    def test_met4d_prepare(self):
        ifq1 = "reads_R1.fastq.gz"
        ifq2 = "reads_R2.fastq.gz"
        ifq1Handle = gzip.open(ifq1, "w")
        ifq2Handle = gzip.open(ifq2, "w")
        
        # pair 1: read 1 has perfect tag and cut site of ind 2 at bp 2; read 2 has the tag but not the cut site
        txt = "@INST1:1:FLOW1:2:2104:15343:197391 1:N:0\n"
        txt += "ATTT" # tag with a 1-bp shift
        txt += "CAGCCCTGGAGTTCCAC\n" # insert with ApekI cut site
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~~\n"
        ifq1Handle.write(txt)
        txt = "@INST1:1:FLOW1:2:2104:15343:197391 2:N:0\n"
        txt += "ATTT" # tag with a 1-bp shift
        txt += "GTAGCTGAGATCGGAAG\n" # insert
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~~\n"
        ifq2Handle.write(txt)
        
        # pair 2: read 1 has perfect tag of ind 1 at bp 2 but wrong cut site; read 2 has tag of ind 1 and the right cut site
        txt = "@INST1:1:FLOW1:2:2104:15343:197392 1:N:0\n"
        txt += "TAAA" # tag with a 1-bp shift
        txt += "TAGCTACATAACTACAT\n" # insert with no ApekI cut site
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~~\n"
        ifq1Handle.write(txt)
        txt = "@INST1:1:FLOW1:2:2104:15343:197392 2:N:0\n"
        txt += "TAAA" # tag with a 1-bp shift
        txt += "CAGCCCTGGAGTTCCAA\n" # insert with ApekI cut site
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~~\n"
        ifq2Handle.write(txt)

        # pair 3: read 1 has perfect tag of ind 3 at bp 2 but wrong cut site; read 2 has wrong tag but right cut site
        txt = "@INST1:1:FLOW1:2:2104:15343:197393 1:N:0\n"
        txt += "TGGG" # tag with a 1-bp shift
        txt += "TAGCTACATTTACTACT\n" # insert with no ApekI cut site
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~~\n"
        ifq1Handle.write(txt)
        txt = "@INST1:1:FLOW1:2:2104:15343:197393 2:N:0\n"
        txt += "TGTG" # tag with a 1-bp shift
        txt += "CAGCCCTGGAGTTCCAG\n" # insert with ApekI cut site
        txt += "+\n"
        txt += "~~~~~~~~~~~~~~~~~~~~~\n"
        ifq2Handle.write(txt)
                
        ifq1Handle.close()
        ifq2Handle.close()
        
        for f in ["test_ind1_R1.fastq.gz", "test_ind1_R2.fastq.gz","test_ind2_R1.fastq.gz", "test_ind2_R2.fastq.gz",
                  "test_unassigned_R1.fastq.gz", "test_unassigned_R2.fastq.gz"]:
            if os.path.isfile(f):
                os.remove(f)
                
        it = "tags.fa"  #Input file with individual ID and tags
        self.writeTagFile(it)
        
        return ifq1, ifq2, it
        
        
    def test_met4d_comp(self, msgs):
        if not os.path.exists("test_ind1_R1.fastq.gz") or \
           not os.path.exists("test_ind1_R2.fastq.gz") or \
           not os.path.exists("test_ind2_R1.fastq.gz") or \
           not os.path.exists("test_ind2_R2.fastq.gz"):
            print("test_met4d: fail (1)")
            return
        else:
            with gzip.open("test_ind1_R1.fastq.gz") as inFqHandle1, \
                 gzip.open("test_ind1_R2.fastq.gz") as inFqHandle2:
                l1 = list(SeqIO.parse(inFqHandle1, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
                l2 = list(SeqIO.parse(inFqHandle2, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
                if len(l1) != 1 or len(l2) != 1:
                    print("test_met4d: fail (2)")
                    return
                if l1[0].id != "INST1:1:FLOW1:2:2104:15343:197392" or \
                   l2[0].id != "INST1:1:FLOW1:2:2104:15343:197392":
                    print("test_met4d: fail (3)")
                    return
                if str(l1[0].seq) != "TAGCTACATAACTACAT" or \
                   str(l2[0].seq) != "CAGCCCTGGAGTTCCAA":
                    print("test_met4d: fail (4)")
                    return
            with gzip.open("test_ind2_R1.fastq.gz") as inFqHandle1, \
                 gzip.open("test_ind2_R2.fastq.gz") as inFqHandle2:
                l1 = list(SeqIO.parse(inFqHandle1, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
                l2 = list(SeqIO.parse(inFqHandle2, "fastq",
                                      alphabet=IUPAC.ambiguous_dna))
                if len(l1) != 1 or len(l2) != 1:
                    print("test_met4d: fail (2)")
                    return
                if l1[0].id != "INST1:1:FLOW1:2:2104:15343:197391" or \
                   l2[0].id != "INST1:1:FLOW1:2:2104:15343:197391":
                    print("test_met4d: fail (3)")
                    return
                if str(l1[0].seq) != "CAGCCCTGGAGTTCCAC" or \
                   str(l2[0].seq) != "GTAGCTGAGATCGGAAG":
                    print("test_met4d: fail (4)")
                    return
        if os.path.exists("test_ind3_R1.fastq.gz") or \
           os.path.exists("test_ind3_R2.fastq.gz"):
            print("test_met4d: fail (5)")
            return
        print("test_met4d: pass")
        
        
    def test_met4d(self):
        if self.verbose > 0:
            print("launch test met4d ...")
            sys.stdout.flush()
        cwd, testDir = self.beforeTest()
        ifq1, ifq2, it = self.test_met4d_prepare()
        msgs = self.launchProg(ifq1, ifq2, it, "4d", 10, "ApeKI", False)
        self.test_met4d_comp(msgs)
        self.afterTest(cwd, testDir)
        
        
    #==========================================================================
    


    
    
    def run(self):
        if "1" in self.testsToRun:
            self.test_met1()
        if "2" in self.testsToRun:
            self.test_met2()
        if "4a" in self.testsToRun:
            self.test_met4a()
        if "4b" in self.testsToRun:
            self.test_met4b()
        if "4c" in self.testsToRun:
            self.test_met4c()
        if "4d" in self.testsToRun:
            self.test_met4d()            
            
if __name__ == "__main__":
    i = TestDemultiplex()
    
    i.setAttributesFromCmdLine()
    
    i.checkAttributes()
    
    if i.verbose > 0:
        startTime = time.time()
        msg = "START %s %s" % (os.path.basename(sys.argv[0]),
                               time.strftime("%Y-%m-%d %H:%M:%S"))
        msg += "\ncmd-line: %s" % ' '.join(sys.argv)
        msg += "\ncwd: %s" % os.getcwd()
        print(msg); sys.stdout.flush()
        
    i.run()
    
    if i.verbose > 0:
        msg = "END %s %s" % (os.path.basename(sys.argv[0]),
                             time.strftime("%Y-%m-%d %H:%M:%S"))
        endTime = time.time()
        runLength = datetime.timedelta(seconds=
                                       math.floor(endTime - startTime))
        msg += " (%s" % str(runLength)
        if "linux" in sys.platform:
            p = Popen(["grep", "VmHWM", "/proc/%s/status" % os.getpid()],
                      shell=False, stdout=PIPE).communicate()
            maxMem = p[0].split()[1]
            msg += "; %s kB)" % maxMem
        else:
            msg += ")"
        print(msg); sys.stdout.flush()
