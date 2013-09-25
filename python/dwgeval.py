import pysam
import re
from bamslider import BamSlider
import sys
from collections import Counter, namedtuple
import pdb


DWGSIM_MATCHER = re.compile(r"(?P<seqname>[\w\-_]+)_(?P<start_1>\d+)_(?P<start_2>\d+)_(?P<strand_1>\d+)_(?P<strand_2>\d+)"
                            "_(?P<random_read_1>\d+)_(?P<random_read_2>\d+)_(?P<nerrors_1>\d+):(?P<nsnps_1>\d+):"
                            "(?P<nindels_1>\d+)_(?P<nerrors_2>\d+):(?P<nsnps_2>\d+):(?P<nindels_2>\d+)_(?P<rnum>\w+)")

DWGSIM_FIELDS = {'seqname':str, 'start_1':int, 'start_2':int, 'strand_1':int, 'strand_2':int, 
                 'random_read_1':int, 'random_read_2':int, 
                 'nerrors_1':int, 'nsnps_1':int, 'nindels_1':int, 
                 'nerrors_2':int, 'nsnps_2':int, 'nindels_2':int, 
                 'rnum':str}

DwgsimRead = namedtuple('DwgsimRead', DWGSIM_FIELDS)

def mean(x):
    return float(sum(x))/len(x) if len(x) > 0 else float('nan')

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

def mapq_slider(window_reads):
    """Take a window (window tuple and reads) and summarize mapq. Return a
    BED entry as a string.
    """
    window, reads = window_reads
    mean_mapq = mean([r.mapq for r in reads])
    mean_mapq = mean_mapq if mean_mapq != float('nan') else "NA"
    bed_line = "\t".join([window.seqname, str(window.start), str(window.end), str(mean_mapq)])
    return bed_line

def aln_slider(window_reads, refs, wiggle=5):
    """Take a window (window tuple and reads) and summarize false positive
    and true positive rates. Return a BED entry as a string.
    """
    assert(wiggle < 20) # just in case a crazy person does something silly
    window, reads = window_reads
    tp = sum([is_correct_aln(r, refs[r.tid], wiggle) for r in reads])    
    bed_line = "\t".join([window.seqname, str(window.start), str(window.end), str(tp), str(len(reads))])
    return bed_line


if __name__ == "__main__":
    size = 1000
    step = 50
    slider = BamSlider(sys.argv[1], size, step)
    refs = slider._samfile.references
    for window_reads in slider:
        print aln_slider(window_reads, refs)

