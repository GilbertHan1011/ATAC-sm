# Tachyon backend pilot: B56-1_ATAC19

Validation root: `/storage/zhangkaiLab/hanlitian/macrophage/validation/20260723_tachyon_backend`.
Tachyon source was `d30dc24`; the glibc-2.34-compatible binary SHA256 was
`896b606e092c2239fac23de4a2061917390da47aa743a11dba73d5664c116e3e`.

The isolated legacy arm disabled chrM prealignment; both arms used the same
declared paired FASTQs. Jobs 2854176 (Tachyon) and the reused legacy outputs
completed on Slurm. Comparison job 2854650 verified both BAMs with
`samtools quickcheck`, coordinate headers, and indexed `chr1:1-100000` queries.

| Metric | Legacy | Tachyon | Result |
|---|---:|---:|---|
| BAM records | 107,212,753 | 107,677,748 | +0.43% |
| chrM fraction | 6.1959% | 6.1964% | +0.0006 pp |
| median fragment | 64 bp | 64 bp | pass |
| fragment histogram Spearman (≤1 kb) | — | 0.99974 | pass |
| peak count | 269,129 | 269,977 | +0.32% |
| reciprocal peak overlap | 97.51% | 97.20% | pass |
| peak bp Jaccard | — | 0.9532 | pass |
| 10 kb BigWig Pearson | — | 0.99989 | pass |
| TSS Pearson | — | 0.9999997 | pass |
| TSS maximum | 4.55009 | 4.54786 | -0.05% |

ATAQV peak, mitochondrial, and TSS metrics differ by less than 1%. Duplicate
rates differ (75.32% legacy, 79.67% Tachyon; +4.34 pp), so the strict duplicate
agreement gate is not met. This is expected to reflect different duplicate
marking implementations and requires expanded-cohort confirmation before any
promotion. The backend remains opt-in.

Machine-readable results: `$validation_root/comparison/atac_pilot.json`.
