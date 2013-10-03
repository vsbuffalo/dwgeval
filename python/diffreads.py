"""
diffreads.py -- return a BAM or BED file of the positions of where the
aligned position of reads don't match the DWGSIM query name. 
"""
import sys
import pysam
import re
from collections import namedtuple

DWGSIM_MATCHER = re.compile(r"(?P<seqname>[\w\-_]+)_(?P<start_1>\d+)_(?P<start_2>\d+)_(?P<strand_1>\d+)_(?P<strand_2>\d+)"
                            "_(?P<random_read_1>\d+)_(?P<random_read_2>\d+)_(?P<nerrors_1>\d+):(?P<nsnps_1>\d+):"
                            "(?P<nindels_1>\d+)_(?P<nerrors_2>\d+):(?P<nsnps_2>\d+):(?P<nindels_2>\d+)_(?P<rnum>\w+)")

DWGSIM_FIELDS = {'seqname':str, 'start_1':int, 'start_2':int, 'strand_1':int, 'strand_2':int, 
                 'random_read_1':int, 'random_read_2':int, 
                 'nerrors_1':int, 'nsnps_1':int, 'nindels_1':int, 
                 'nerrors_2':int, 'nsnps_2':int, 'nindels_2':int, 
                 'rnum':str}

DwgsimRead = namedtuple('DwgsimRead', DWGSIM_FIELDS)

def dwgsim_parser(qname):
    """Parse the header for dwgsim simulated reads.
    """
    matches = DWGSIM_MATCHER.match(qname)
    assert(matches is not None)
    groups = dict([(k, DWGSIM_FIELDS[k](v)) for k, v in matches.groupdict().items()])
    return DwgsimRead(**groups)

def is_correct_aln(read, rname, wiggle=5):
    """See if a read is aligned within `wiggle` of it's real start base.

    Note: we just compare the read.pos to *both* positions in the
    header (TODO, maybe be a bit more explicit). This will be fine
    though because our fragment sizes will never be small enough (less
    than wiggle!) to screw this up.
    """
    dwgsim_read = dwgsim_parser(read.qname)
    if rname != dwgsim_read.seqname:
        return False
    return (dwgsim_read.start_1 - read.pos) <= wiggle or (dwgsim_read.start_2 - read.pos) <= wiggle    

if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit("usage: bam/bed wiggle file.bam")

    format = sys.argv[1]
    wiggle = int(sys.argv[2])
    samfile = pysam.Samfile(sys.argv[3])
    refs = samfile.references

    if format == "bam":
        outbamfile = pysam.Samfile("-", "wb", template=samfile)
    elif format == "bed":
        pass
    else:
        sys.exit("choose bed or bam")
        
    total, incorrect, unmapped = 0, 0, 0
    for read in samfile:
        total += 1
        if read.is_unmapped:
            unmapped += 1
            continue
        if not is_correct_aln(read, refs[read.tid], wiggle=wiggle):
            if format == "bam":
                outbamfile.write(read)
            else:
                dread = dwgsim_parser(read.qname)
                assert(dread.strand_1 != dread.strand_2)
                fstrand = dread.strand_1 == 0
                # assumes pairs have each have same length!
                rlen = len(read.seq)
                start = dread.start_1 if fstrand else dread.start_2
                end = dread.start_2+rlen if fstrand else dread.start_1+rlen
                assert(start < end)
                fields = [dread.seqname, start, end]
                print "\t".join(map(str, fields))
    sys.stderr.write("incorrect: %d\nunmapped: %d\ntotal: %d\n" % (incorrect, unmapped, total))
    
        
            
        

