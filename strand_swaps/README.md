# Install

```shell
cd strand_swaps/
pip install -e .
```

# Use

```shell
strand_swaps run --input file.gbk
```

# Output

You get histograms:
```text
file.consec_strand.tsv - counts of consecutive genes on same strand
file.consec_frame.tsv - counts of consecutive genes on same frame
file.consec_overlap.tsv - counts of consecutive genes that overlap
```
