#!/usr/bin/env python

import os
import sys
import pysam

class Barcode(object):
    outfile = ""
    out = None
    count = 0

    def __init__(self, outfile, template):
        self.outfile = outfile
        self.out = pysam.AlignmentFile(outfile, "wb", template=template)
        self.count = 0

    def write(self, read):
        self.count += 1
        self.out.write(read)

    def close(self):
        self.out.close()

class Main(object):
    barcodesfile = None
    bamfile = None
    infile = None
    bctable = {}
    prefix = ""
    reportfile = None
    delEmpty = False
    
    def usage(self):
        sys.stdout.write("""splitbam.py - write reads to different BAM files on the basis of their RG tag.

Usage: splitbam.py [options] barcodes bamfile

This program reads a list of barcodes from the `barcodes' file (one per line), then reads the bamfile and
checks the RG tag of each read. If the RG tag is in the provided list of barcodes, the read is written to
a new bamfile called `prefix-barcode.bam', where `prefix' can be specified with the -p option.  A final
report listing the number of reads associated with each barcode can be generated with the -r option.

Options:

  -p P : Use P as the prefix for output files (default: {}).
  -r R : Write final report to file R.
  -d   : Delete empty BAM files when done.

""".format(self.prefix))

    def parseArgs(self, args):
        if "-h" in args or "--help" in args:
            return False
        prev = ""
        for a in args:
            if prev == "-p":
                self.prefix = a
                prev = ""
            elif prev == "-r":
                self.reportfile = a
                prev = ""
            elif a in ["-p", "-r"]:
                prev = a
            elif a == "-d":
                self.delEmpty = True
            elif self.barcodesfile is None:
                self.barcodesfile = a
            else:
                self.bamfile = a
        return self.barcodesfile and self.bamfile

    def readBarcodes(self):
        self.bctable = {}
        with open(self.barcodesfile, "r") as f:
            for line in f:
                bc = line.strip()
                #outfile = "{}-{}.bam".format(self.prefix, bc)
                outfile = "{}.bam".format(bc)
                self.bctable[bc] = Barcode(outfile, self.infile)

    def run(self):
        sys.stderr.write("*** Splitting BAM file {}\n".format(self.bamfile))
        self.infile = pysam.AlignmentFile(self.bamfile, "rb")
        self.readBarcodes()
        sys.stderr.write("*** {} barcodes read from {}\n".format(len(self.bctable), self.barcodesfile))
        nin = nout = 0
        for read in self.infile:
            nin += 1
            rg = read.get_tag("RG")
            if rg in self.bctable:
                self.bctable[rg].write(read)
                nout += 1
#            if nin % 1000 == 0:
#                sys.stderr.write("\r{}".format(nin))
        sys.stderr.write("*** Closing all BAM files.\n")
        for bc in self.bctable:
            bam = self.bctable[bc]
            bam.close()
            if self.delEmpty and bam.count == 0:
                os.remove(bam.outfile)
                
        if self.reportfile:
            with open(self.reportfile, "w") as out:
                for bc in self.bctable:
                    out.write("{}\t{}\n".format(bc, self.bctable[bc].count))
        sys.stderr.write("{} / {} reads written ({:.1f}%)\n".format(nout, nin, 100.0 * nout / nin))

def main(args):
    M = Main()
    if M.parseArgs(args):
        M.run()
    else:
        M.usage()

if __name__ == "__main__":
    main(sys.argv[1:])
