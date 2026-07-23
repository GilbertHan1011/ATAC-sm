
aligner = config["alignment"].get("tool", "bowtie2")
_ALIGNER_LOG_SUFFIX = "txt" if aligner == "bowtie2" else "bwa.log"

# Check if prealignments are enabled
prealign_enabled = config["alignment"].get("prealign", {}).get("enabled", True)
prealignments = config["alignment"].get("prealign", {}).get("indices", []) or []
has_prealignments = prealign_enabled and len(prealignments) > 0

# Resolved after processing.smk has run; defaults to legacy if processing.smk
# did not define it (e.g., a partial include).
ALIGNMENT_BACKEND = config["alignment"].get("backend", "legacy")

# When Tachyon produces the canonical BAM, the legacy producer inputs are not
# available, so the legacy-only MultiQC inputs are dropped instead of pointing
# at fabricated files.
_legacy_only_inputs = (
    ALIGNMENT_BACKEND == "tachyon_upstream"
)
_fastp_html = [] if _legacy_only_inputs else expand(
    os.path.join(result_path, 'report', 'fastp', '{sample_run}_fastp.html'),
    sample_run=annot.index.tolist(),
)
_fastp_json = [] if _legacy_only_inputs else expand(
    os.path.join(result_path, 'report', 'fastp', '{sample_run}_fastp.json'),
    sample_run=annot.index.tolist(),
)
_aligner_log = [] if _legacy_only_inputs else expand(
    os.path.join(result_path, 'report', 'align', '{sample}.' + _ALIGNER_LOG_SUFFIX),
    sample=samples.keys(),
)
_samblaster_log = [] if _legacy_only_inputs else expand(
    os.path.join(result_path, 'report', 'align', '{sample}.samblaster.log'),
    sample=samples.keys(),
)
_flagstat_log = [] if _legacy_only_inputs else expand(
    os.path.join(result_path, 'report', 'align', '{sample}.samtools_flagstat.log'),
    sample=samples.keys(),
)

rule multiqc:
    input:
        # Legacy-only inputs (fastp reports + legacy aligner logs). Under
        # alignment.backend='tachyon_upstream' these expand to [] so the DAG
        # closes without depending on missing files.
        _fastp_html,
        _fastp_json,
        # Collect per-sample stats and logs for MultiQC
        expand(os.path.join(result_path, 'report', 'align_stats', '{sample}.align.stats.tsv'), sample=samples.keys()),
        _aligner_log,
        _samblaster_log,
        _flagstat_log,
        expand(os.path.join(result_path, 'report', 'tss_coverage', '{sample}.tss_histogram.csv'), sample=samples.keys()),
        sample_annotation = annotation_sheet_path,
        # Collect prealign stats if enabled (per sample, since prealignment now happens at sample level)
        prealign_stats = expand(os.path.join(result_path, 'report',"prealigned", '{sample}.prealign.stats.tsv'), sample=samples.keys()) if has_prealignments else [],
    output:
        multiqc_report = report(os.path.join(result_path,"report","multiqc_report.html"),
                                caption="../report/multiqc.rst",
                                category="{}_{}".format(config["project"]["name"], module_name),
                                subcategory="QC",
                                labels={
                                    "name": "MultiQC report",
                                    "type": "HTML",
                                }),
        multiqc_stats = os.path.join(result_path, "report", "multiqc_report_data", "multiqc_general_stats.txt"),
    params:
        result_path = result_path,
        multiqc_configs = "{{'title': '{name}', 'intro_text': 'Quality Control Metrics of the ATAC-seq pipeline.', 'annotation': '{annot}', 'genome': '{genome}', 'exploratory_columns': {exploratory_columns}, 'skip_versions_section': true}}".format(name = config["project"]["name"], annot = annotation_sheet_path, genome = config["project"]["genome"], exploratory_columns = config["project"].get("annot_columns", "[]"))
    resources:
        mem_mb=config["resources"].get("mem_mb", 16000),
        runtime = 20,
    threads: config["resources"].get("threads", 2)
    conda:
        "../envs/multiqc.yaml",
    log:
        "logs/rules/multiqc.log"
    shell:
        """
        multiqc {params.result_path}/report --force --verbose --outdir {params.result_path}/report --filename multiqc_report.html
        """
