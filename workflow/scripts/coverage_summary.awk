# Summarise `samtools depth -a` output (chrom, pos, depth).
# Variables (pass with -v): sample, low (depth threshold), bed (output path)
# stdout: one-row TSV summary; `bed`: runs of positions with depth < low.
BEGIN {
    OFS = "\t"
    printf "" > bed            # make sure the BED exists even if empty
}
{
    n++
    sum += $3
    if ($3 >= 1)  c1++
    if ($3 >= 10) c10++
    if ($3 >= 30) c30++
    if ($3 < low) {
        if (inrun && $1 == chr && $2 == end + 1) {
            end = $2
        } else {
            if (inrun) print chr, start - 1, end > bed
            chr = $1; start = $2; end = $2; inrun = 1
        }
    } else if (inrun) {
        print chr, start - 1, end > bed
        inrun = 0
    }
}
END {
    if (inrun) print chr, start - 1, end > bed
    print "sample", "genome_positions", "mean_depth", "pct_genome_ge1x", "pct_genome_ge10x", "pct_genome_ge30x"
    if (n == 0) { print sample, 0, 0, 0, 0, 0; exit }
    printf "%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n", sample, n, sum / n, 100 * c1 / n, 100 * c10 / n, 100 * c30 / n
}
