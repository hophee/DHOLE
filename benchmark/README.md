# Benchmark DHOLE

Run from the repository root on Linux/WSL with the DHOLE environment installed:

```bash
bash benchmark/run.sh             # N=10 runs, M=10 CDS, seed=1
bash benchmark/run.sh 20 15 42     # N M seed
```

Requires `test/MG1655.fna`, `test/MG1655.gff`, Bash, tar, and the ggplot2 package
in the `oligo_design` environment. Uses `tools/run-r`; dependencies are not
installed automatically. Synthetic pTarget/pCas fixtures match `test/test_run.sh`.

Each run receives M distinct randomly selected CDS from the GFF (locus_tag,
otherwise gene). CDS may recur across runs. The seed and selected CDS are saved.
DHOLE runs sequentially in separate processes with `--filtering-level 2`;
all other design parameters retain their defaults, including the thread count.
The arm size limits from the integration test are not applied.

Results are created in `benchmark/results-DATE-TIME-PID/`:

- `analysis.pdf`: ggplot charts of runtime (boxplots, individual points, means),
  successful/unsuccessful design percentages, Tm and GC histograms, and bar plots
  of candidate filtering categories, rejection reasons, and failure stages.
- `timings.tsv`, `mean_times.tsv`: per-run timings and arithmetic means.
  Real is the elapsed time of the entire DHOLE process; CPU is user + system time,
  including child processes that are waited for. Initialization and index building
  are included; CDS sampling, analysis, and archiving are excluded. Failed runs
  are included.
- `targets.tsv`, `success_rates.tsv`: success means `status=ok` in
  `design_summary.tsv`; a missing row also counts as a failure. The denominator
  is all N×M attempts. Process exit codes are recorded separately: a zero exit
  code does not imply a successful CDS design.
- `primers.tsv`: the six final PCR primers from each successful design.
  Tm comes from the DHOLE WetLab report (rounded to 0.1 °C); GC is calculated
  from the annealing sequence without auxiliary tails. sgRNA oligonucleotides
  are excluded. Repeated CDS are counted as separate attempts.
- `rankings.tsv`, `filter_counts.tsv`, `rejection_counts.tsv`: all available
  candidate rows, including failed targets; stages and attempts may count the
  same pair more than once. Unevaluated filters are listed separately;
  multiple rejection reasons may apply to one candidate.
- `sampled_cds.tsv`, `metadata.txt`: samples, seed, input checksums, and R/package
  versions; DHOLE parameters are stored in its `run_parameters.tsv`.
- `dhole-results.tar.gz`: created automatically after analysis; contains `raw/`
  with all DHOLE outputs, plasmids, commands, and run logs. The original `raw/`
  directory is preserved; the archive and PDF are saved alongside it.

If no data are available for a plot, the PDF includes a corresponding note.
A DHOLE crash does not stop subsequent runs. Analysis or archiving errors cause
the benchmark to exit with an error while preserving the raw data.

Test the analysis without running DHOLE:

```bash
bash tools/run-r benchmark/test_analysis.R
```
