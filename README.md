# dwgeval -- Tools to evaluate DWGSIM alignements.

[dwgsim](https://github.com/nh13/DWGSIM) is a read simulation
tool. This set of tools makes working with BAM/SAM files that are
aligned dwgsim reads easier. Dwgsim gives the read position of the
read in the query name.

Note that there are two things we may want to see (which are not
included in dwgsim's `dwgsim_eval` output). Correct alignment will
always be the difference between the real position in the query name
and the alignment position, but where the read is placed will affect
what we tell using downstream analysis:

1. Where misaligned reads are aligning to.

2. Where misaligned reads are coming from.

For high-coverage simulated data, these can be tricky to calculate
quickly. There are two approaches: for (1) for a sorted BAM file, use
sliding windows to calculate misalignment statistics across alignment
windows. For (2), make this easier by taking all misaligned reads and
converting their real positions (from the query name) into a BED
file.

Note that for windowed-processing, we rely on sorted values a
lot. This is easy for (1) as samtools can easily sort a BAM file, but
sorting by positions in the query name for (2) is less easy, so we
resort to using `python/diffreads.py` to change incorrectly aligned
reads to a BED file, which can then be processed downstread with
`bedtools` or another tool. For (1) we use
[https://github.com/vsbuffalo/bamslider](bamslider) and
`python/dwgeval.py`.


