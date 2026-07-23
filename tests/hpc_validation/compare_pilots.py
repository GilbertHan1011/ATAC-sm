#!/usr/bin/env python3
import csv
import json
import math
import statistics
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pyBigWig


ROOT = Path(sys.argv[1])
SAMTOOLS = sys.argv[2]
OUT = ROOT / "comparison"


def command(*args: str) -> int:
    return int(subprocess.check_output(args, text=True).strip())


def bam_metrics(bam: Path) -> dict[str, float | int | str]:
    subprocess.run([SAMTOOLS, "quickcheck", "-v", str(bam)], check=True)
    header = subprocess.check_output([SAMTOOLS, "view", "-H", str(bam)], text=True)
    histogram: Counter[int] = Counter()
    with subprocess.Popen(
        [SAMTOOLS, "view", "-f", "66", str(bam)], text=True, stdout=subprocess.PIPE
    ) as proc:
        assert proc.stdout
        for line in proc.stdout:
            fields = line.split("\t")
            if len(fields) > 8:
                length = abs(int(fields[8]))
                if 0 < length <= 1000:
                    histogram[length] += 1
        if proc.wait():
            raise RuntimeError(f"fragment scan failed: {bam}")
    midpoint = (sum(histogram.values()) - 1) // 2
    seen = 0
    median = 0
    for length, count in sorted(histogram.items()):
        seen += count
        if seen > midpoint:
            median = length
            break
    return {
        "quickcheck": True,
        "sort_order": next((x for x in header.splitlines() if x.startswith("@HD")), ""),
        "records": command(SAMTOOLS, "view", "-c", str(bam)),
        "chr1_region_records": command(SAMTOOLS, "view", "-c", str(bam), "chr1:1-100000"),
        "unmapped_records": command(SAMTOOLS, "view", "-c", "-f", "4", str(bam)),
        "duplicate_records": command(SAMTOOLS, "view", "-c", "-f", "1024", str(bam)),
        "chrM_records": command(SAMTOOLS, "view", "-c", str(bam), "chrM"),
        "median_fragment_length": median,
        "fragment_histogram": dict(histogram),
    }


def pearson(left: list[float], right: list[float]) -> float:
    pairs = [(x, y) for x, y in zip(left, right) if not (math.isnan(x) or math.isnan(y))]
    if len(pairs) < 2:
        return float("nan")
    xs, ys = zip(*pairs)
    xmean, ymean = statistics.mean(xs), statistics.mean(ys)
    numerator = sum((x - xmean) * (y - ymean) for x, y in pairs)
    denom = math.sqrt(sum((x - xmean) ** 2 for x in xs) * sum((y - ymean) ** 2 for y in ys))
    return numerator / denom if denom else float("nan")


def rank(values: list[float]) -> list[float]:
    ordered = sorted(enumerate(values), key=lambda pair: pair[1])
    result = [0.0] * len(values)
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and ordered[end][1] == ordered[start][1]:
            end += 1
        value = (start + end - 1) / 2
        for index, _ in ordered[start:end]:
            result[index] = value
        start = end
    return result


def histogram_spearman(left: dict[str, int], right: dict[str, int]) -> float:
    keys = sorted({int(key) for key in left} | {int(key) for key in right})
    return pearson(
        rank([left.get(str(key), left.get(key, 0)) for key in keys]),
        rank([right.get(str(key), right.get(key, 0)) for key in keys]),
    )


def read_peaks(path: Path) -> list[tuple[str, int, int]]:
    peaks = []
    with path.open() as handle:
        for line in handle:
            chrom, start, end, *_ = line.split("\t")
            peaks.append((chrom, int(start), int(end)))
    return sorted(peaks)


def interval_overlap(left: list[tuple[str, int, int]], right: list[tuple[str, int, int]]) -> dict[str, float]:
    left_hit, right_hit, shared_bp = set(), set(), 0
    i = j = 0
    while i < len(left) and j < len(right):
        lc, ls, le = left[i]
        rc, rs, re = right[j]
        if lc < rc or (lc == rc and le <= rs):
            i += 1
        elif rc < lc or (lc == rc and re <= ls):
            j += 1
        else:
            shared_bp += max(0, min(le, re) - max(ls, rs))
            left_hit.add(i)
            right_hit.add(j)
            if le <= re:
                i += 1
            else:
                j += 1
    left_bp = sum(end - start for _, start, end in left)
    right_bp = sum(end - start for _, start, end in right)
    return {
        "legacy_peak_count": len(left),
        "tachyon_peak_count": len(right),
        "legacy_reciprocal_overlap": len(left_hit) / len(left) if left else 0,
        "tachyon_reciprocal_overlap": len(right_hit) / len(right) if right else 0,
        "peak_bp_jaccard": shared_bp / (left_bp + right_bp - shared_bp) if left_bp + right_bp > shared_bp else 0,
    }


def bigwig_pearson(left: Path, right: Path) -> float:
    with pyBigWig.open(str(left)) as first, pyBigWig.open(str(right)) as second:
        values_left, values_right = [], []
        for chrom, length in first.chroms().items():
            if chrom not in second.chroms():
                continue
            for start in range(0, min(length, second.chroms()[chrom]), 10_000):
                end = min(start + 10_000, length, second.chroms()[chrom])
                values_left.append(first.stats(chrom, start, end, type="mean")[0] or 0.0)
                values_right.append(second.stats(chrom, start, end, type="mean")[0] or 0.0)
    return pearson(values_left, values_right)


def tsv_percent(path: Path) -> float:
    with path.open() as handle:
        return float(next(csv.DictReader(handle, delimiter="\t"))["percent"])


def ataqv(path: Path) -> dict[str, float | int]:
    metrics = json.loads(path.read_text())[0]["metrics"]
    keys = ["paired_reads", "duplicate_reads", "total_mitochondrial_reads", "ppm_in_peaks", "total_peaks", "tss_enrichment"]
    return {key: metrics[key] for key in keys}


def tss_pearson(left: Path, right: Path) -> dict[str, float]:
    def values(path: Path) -> list[float]:
        with path.open() as handle:
            return [float(row["count"]) for row in csv.DictReader(handle)]
    a, b = values(left), values(right)
    return {"tss_pearson": pearson(a, b), "legacy_tss_max": max(a), "tachyon_tss_max": max(b)}


def write_json(name: str, data: dict) -> None:
    (OUT / f"{name}.json").write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def compare_atac() -> None:
    base = ROOT / "results" / "atac"
    sample = "B56-1_ATAC19"
    legacy, tachyon = base / "legacy_no_prealign", base / "tachyon"
    data = {
        "legacy": bam_metrics(legacy / "important_processed/bam" / f"{sample}.filtered.bam"),
        "tachyon": bam_metrics(tachyon / "important_processed/bam" / f"{sample}.filtered.bam"),
    }
    data["fragment_histogram_spearman"] = histogram_spearman(data["legacy"]["fragment_histogram"], data["tachyon"]["fragment_histogram"])
    data["peaks"] = interval_overlap(
        read_peaks(legacy / "important_processed/peaks" / f"{sample}_peaks.narrowPeak"),
        read_peaks(tachyon / "important_processed/peaks" / f"{sample}_peaks.narrowPeak"),
    )
    data["bigwig_10kb_pearson"] = bigwig_pearson(legacy / "important_processed/tracks" / f"{sample}.bw", tachyon / "important_processed/tracks" / f"{sample}.bw")
    data["tss"] = tss_pearson(legacy / "report/tss_coverage" / f"{sample}.tss_histogram.csv", tachyon / "report/tss_coverage" / f"{sample}.tss_histogram.csv")
    data["ataqv"] = {"legacy": ataqv(legacy / "report/ataqv" / f"{sample}.ataqv.json"), "tachyon": ataqv(tachyon / "report/ataqv" / f"{sample}.ataqv.json")}
    data["tachyon_stats"] = json.loads((tachyon / "report/align" / f"{sample}.tachyon.stats.json").read_text())
    write_json("atac_pilot", data)


def compare_cuttag() -> None:
    base = ROOT / "results" / "cuttag"
    sample = "B51-1_CUT22"
    legacy, tachyon = base / "legacy", base / "tachyon"
    data = {
        "legacy": bam_metrics(legacy / "Important_processed/Bam" / f"{sample}.sorted.markd.bam"),
        "tachyon": bam_metrics(tachyon / "Important_processed/Bam" / f"{sample}.sorted.markd.bam"),
    }
    data["fragment_histogram_spearman"] = histogram_spearman(data["legacy"]["fragment_histogram"], data["tachyon"]["fragment_histogram"])
    data["peaks"] = interval_overlap(
        read_peaks(legacy / "Important_processed/Peaks/callpeaks" / f"macs2_narrow_{sample}_peaks.narrowPeak"),
        read_peaks(tachyon / "Important_processed/Peaks/callpeaks" / f"macs2_narrow_{sample}_peaks.narrowPeak"),
    )
    data["bigwig_10kb_pearson"] = bigwig_pearson(legacy / "Important_processed/Track/tracks" / f"{sample}.bw", tachyon / "Important_processed/Track/tracks" / f"{sample}.bw")
    data["frip_percent"] = {"legacy": tsv_percent(legacy / "Report/plotEnrichment" / f"frip_{sample}.tsv"), "tachyon": tsv_percent(tachyon / "Report/plotEnrichment" / f"frip_{sample}.tsv")}
    data["tachyon_stats"] = json.loads((tachyon / "Report/tachyon_upstream" / f"{sample}.stats.json").read_text())
    write_json("cuttag_pilot", data)


OUT.mkdir(exist_ok=True)
compare_atac()
compare_cuttag()
