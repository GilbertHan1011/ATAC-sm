# ============================================================================
# Conditional aligner selection using if-else
# Only ONE aligner rule will be defined based on config
# ============================================================================

if config.get("aligner", "bowtie2") == "bowtie2":
    # alignment with bowtie2 & samtools (per run)
    rule align_bowtie2:
        input:
            fasta_fwd = os.path.join(result_path, "trimmed", "{sample_run}_1.fq.gz"),
            fasta_rev = os.path.join(result_path, "trimmed", "{sample_run}_2.fq.gz"),
            bowtie2_index = os.path.dirname(config["bowtie2_index"]),
            adapter_fasta = config["adapter_fasta"] if config["adapter_fasta"]!="" else [],
            whitelisted_regions = config["whitelisted_regions"],
        wildcard_constraints:
            sample_run="|".join(annot.index.tolist())  # Only match actual sample_run names from annotation
        output:
            bam = temp(os.path.join(result_path,"results","{sample_run}","mapped", "{sample_run}.bam")),
            output_bai =  temp(os.path.join(result_path,"results","{sample_run}","mapped", "{sample_run}.bam.bai")),
            filtered_bam = os.path.join(result_path,"results","{sample_run}","mapped", "{sample_run}.filtered.bam"),
            filtered_bai = os.path.join(result_path,"results","{sample_run}","mapped", "{sample_run}.filtered.bam.bai"),
            bowtie_log = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.txt'),
            bowtie_met = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.bowtie2.met'),
            samblaster_log = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.samblaster.log'),
            samtools_log = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.samtools.log'),
            samtools_flagstat_log = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.samtools_flagstat.log'),
            stats = os.path.join(result_path, 'results', "{sample_run}", '{sample_run}.align.stats.tsv'),
        params:
            sample_name = lambda w: annot.loc[w.sample_run, 'sample_name'],
            bowtie2_input = lambda w, input: f"-1 {input.fasta_fwd} -2 {input.fasta_rev}" if annot.loc[w.sample_run, "read_type"] == "paired" else f"-U {input.fasta_fwd}",
            filtering = lambda w: "-q 30 -F {flag} -f 2 -L {whitelist}".format(flag=config['SAM_flag'], whitelist=config["whitelisted_regions"]) if annot.loc[w.sample_run, "read_type"] == "paired" else "-q 30 -F {flag} -L {whitelist}".format(flag=config['SAM_flag'], whitelist=config["whitelisted_regions"]),
            add_mate_tags = lambda w: "--addMateTags" if annot.loc[w.sample_run, "read_type"] == "paired" else " ",
            adapter_sequence = "-a " + config["adapter_sequence"] if config["adapter_sequence"] !="" else " ",
            adapter_fasta = "--adapter_fasta " + config["adapter_fasta"] if config["adapter_fasta"] !="" else " ",
            sequencing_platform = config["sequencing_platform"],
            mitochondria_name = config["mitochondria_name"],
            bowtie2_index = config["bowtie2_index"], # The basename of the index for the reference genome excluding the file endings e.g., *.1.bt2
            bowtie2_local_mode = "--local" if config.get("local", False) else "",  # Add this param for local mode
        resources:
            mem_mb=config.get("mem", "16000"),
        threads: 4*config.get("threads", 2)
        conda:
            "../envs/bowtie2.yaml",
        log:
            "logs/rules/align_{sample_run}.log"
        shell:
            """
            mkdir -p $(dirname {output.stats})
            result_path=$(dirname {output.stats})
            find $result_path -type f -name '*.bam.tmp.*' -exec rm {{}} +;
            
            RG="--rg-id {wildcards.sample_run} --rg SM:{params.sample_name} --rg PL:{params.sequencing_platform}"

            bowtie2 $RG --very-sensitive --no-discordant -p {threads} --maxins 2000 -x {params.bowtie2_index} \
                {params.bowtie2_local_mode}  \
                --met-file "{output.bowtie_met}" {params.bowtie2_input} 2> "{output.bowtie_log}" | \
                samblaster {params.add_mate_tags} 2> "{output.samblaster_log}" | \
                samtools sort -o "{output.bam}" - 2>> "{output.samtools_log}";
                
            
            samtools index "{output.bam}" 2>> "{output.samtools_log}";
            samtools idxstats "{output.bam}" | awk '{{ sum += $3 + $4; if($1 == "{params.mitochondria_name}") {{ mito_count = $3; }}}}END{{ print "mitochondrial_fraction\t"mito_count/sum }}' > "{output.stats}";
            samtools flagstat "{output.bam}" > "{output.samtools_flagstat_log}";

            samtools view {params.filtering} -o "{output.filtered_bam}" "{output.bam}";
            samtools index "{output.filtered_bam}";
          
            """

elif config.get("aligner", "bowtie2") == "bwa-mem2":
    def _sanitize_bwa_min_score_flag():
        value = config.get("bwa_min_score")
        if value is None:
            return ""
        if isinstance(value, str):
            stripped = value.strip()
            if stripped.lower() in {"", "null", "none"}:
                return ""
            value_str = stripped
        else:
            value_str = value
        return f"-T {value_str}"

    def _resolve_bwa_index_prefix():
        """
        Resolve the path prefix to pass to bwa/bwa-mem2.
        If config['bwa_mem2_path'] points directly to a prefix (i.e., *.0123 exists),
        use it as-is. If it points to a directory, try to join with the genome_fasta
        stem (basename without extension).
        """
        import os
        idx = config.get("bwa_mem2_path")
        if not idx or idx in ("", "null"):
            # fallback to genome_fasta (bwa classic index prefix)
            return config["genome_fasta"]
        # If directly a prefix (files exist), use as-is
        if any(os.path.exists(f"{idx}{ext}") for ext in (".0123", ".amb", ".ann", ".bwt", ".pac", ".sa")):
            return idx
        # If a directory, try to derive prefix from genome_fasta stem
        if os.path.isdir(idx):
            stem = os.path.splitext(os.path.basename(config["genome_fasta"]))[0]
            candidate = os.path.join(idx, stem)
            if any(os.path.exists(f"{candidate}{ext}") for ext in (".0123", ".amb", ".ann", ".bwt", ".pac", ".sa")):
                return candidate
        # As a last resort, return provided value (will likely fail fast and inform the user)
        return idx

    # alignment with BWA MEM & samtools (per run)
    rule align_bwa_mem:
        input:
            fasta_fwd = lambda w: os.path.join(result_path, "trimmed", f"{w.sample_run}_1.fq.gz") if annot.loc[w.sample_run, "read_type"] == "paired" else os.path.join(result_path, "trimmed", f"{w.sample_run}.fq.gz"),
            fasta_rev = lambda w: os.path.join(result_path, "trimmed", f"{w.sample_run}_2.fq.gz") if annot.loc[w.sample_run, "read_type"] == "paired" else [],
            # If bwa_mem2_path is set, use BWA-MEM2 index directory (pre-indexed)
            # If bwa_mem2_path is not set, require bwa_index rule to create index files
            index = lambda w: multiext(config["genome_fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa", ".0123", ".alt") if not config.get("bwa_mem2_path") or config.get("bwa_mem2_path") == "null" or config.get("bwa_mem2_path") == "" else [],
            whitelisted_regions = config["whitelisted_regions"],
        wildcard_constraints:
            sample_run="|".join(annot.index.tolist())  # Only match actual sample_run names from annotation
        output:
            bam = temp(os.path.join(result_path,"results","{sample_run}","mapped", "{sample_run}.bam")),
            output_bai =  temp(os.path.join(result_path,"results","{sample_run}","mapped", "{sample_run}.bam.bai")),
            filtered_bam = os.path.join(result_path,"results","{sample_run}","mapped", "{sample_run}.filtered.bam"),
            filtered_bai = os.path.join(result_path,"results","{sample_run}","mapped", "{sample_run}.filtered.bam.bai"),
            bwa_log = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.bwa.log'),
            samblaster_log = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.samblaster.log'),
            samtools_log = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.samtools.log'),
            samtools_flagstat_log = os.path.join(result_path, 'results', "{sample_run}", 'mapped', '{sample_run}.samtools_flagstat.log'),
            stats = os.path.join(result_path, 'results', "{sample_run}", '{sample_run}.align.stats.tsv')
        params:
            sample_name = lambda w: annot.loc[w.sample_run, 'sample_name'],
            bwa_input = lambda w, input: f"{input.fasta_fwd} {input.fasta_rev}" if annot.loc[w.sample_run, "read_type"] == "paired" else f"{input.fasta_fwd}",
            filtering = lambda w: "-q 30 -F {flag} -f 2 -L {whitelist}".format(flag=config['SAM_flag'], whitelist=config["whitelisted_regions"]) if annot.loc[w.sample_run, "read_type"] == "paired" else "-q 30 -F {flag} -L {whitelist}".format(flag=config['SAM_flag'], whitelist=config["whitelisted_regions"]),
            add_mate_tags = lambda w: "--addMateTags" if annot.loc[w.sample_run, "read_type"] == "paired" else " ",
            sequencing_platform = config["sequencing_platform"],
            mitochondria_name = config["mitochondria_name"],
            bwa_args = config.get("bwa_args", ""),
            bwa_min_score_flag = _sanitize_bwa_min_score_flag(),
            bwa_m_flag = "-M",  # Mark shorter split hits as secondary
            bwa_mem2_index_path = config.get("bwa_mem2_path"),
            # Use bwa-mem2 index path if set, otherwise use genome_fasta (for bwa)
            bwa_index_path = lambda w: _resolve_bwa_index_prefix(),
        resources:
            mem_mb=config.get("mem", "64000"),
        threads: 4*config.get("threads", 2)
        conda:
            "../envs/bwa.yaml",
        log:
            "logs/rules/align_bwa_mem_{sample_run}.log"
        shell:
            """
            mkdir -p $(dirname {output.stats})
            result_path=$(dirname {output.stats})
            find $result_path -type f -name '*.bam.tmp.*' -exec rm {{}} +;
            
            RG="@RG\\tID:{wildcards.sample_run}\\tSM:{params.sample_name}\\tPL:{params.sequencing_platform}"

            BWA_INDEX_PATH="{params.bwa_index_path}"
            
            # Determine which BWA binary to use based on bwa_mem2_path
            if [ -n "{params.bwa_mem2_index_path}" ] && [ "{params.bwa_mem2_index_path}" != "null" ] && [ "{params.bwa_mem2_index_path}" != "" ]; then
                BWA_BINARY="bwa-mem2"
            else
                BWA_BINARY="bwa"
            fi

            $BWA_BINARY mem \
                {params.bwa_args} \
                {params.bwa_m_flag} \
                {params.bwa_min_score_flag} \
                -R "$RG" \
                -t {threads} \
                $BWA_INDEX_PATH \
                {params.bwa_input} 2> {output.bwa_log} | \
                samblaster {params.add_mate_tags} 2> {output.samblaster_log} | \
                samtools view -bhS -F 0x0100 -O BAM - 2>> {output.samtools_log} | \
                samtools sort -o {output.bam} - 2>> {output.samtools_log};
                
            samtools index "{output.bam}" 2>> "{output.samtools_log}";
            samtools idxstats "{output.bam}" | awk '{{ sum += $3 + $4; if($1 == "{params.mitochondria_name}") {{ mito_count = $3; }}}}END{{ print "mitochondrial_fraction\\t"mito_count/sum }}' > "{output.stats}";
            samtools flagstat "{output.bam}" > "{output.samtools_flagstat_log}";

            samtools view {params.filtering} -o "{output.filtered_bam}" "{output.bam}";
            samtools index "{output.filtered_bam}";
            """

else:
    raise ValueError(f"Unknown aligner: {config.get('aligner', 'bowtie2')}. Must be 'bowtie2' or 'bwa-mem2'")

rule bwa_index:
    input:
        fasta = config["genome_fasta"]
    output:
        index = multiext(config["genome_fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa", ".0123", ".alt"),
    log:
        "logs/bwa_index/bwa_index.log"
    conda:
        "../envs/bwa_mem.yaml"
    shell:
        """
        # Create BWA index
        bwa index {input.fasta} 2> {log}
        """

# Merge BAM files for samples with multiple runs
# Note: This rule only matches actual sample names (not sample_run names)
def get_merge_bam_inputs(wildcards):
    """Get inputs for merge_bam, validating that wildcards.sample is a valid sample name"""
    if wildcards.sample not in samples:
        # Return a non-existent file to prevent rule matching for invalid sample names
        # This file will never exist, so the rule won't match
        return [os.path.join(result_path, "__INVALID_SAMPLE__", f"{wildcards.sample}.bam")]
    return [os.path.join(result_path, "results", sr, "mapped", f"{sr}.filtered.bam") 
            for sr in get_runs_for_sample(wildcards.sample)]

rule merge_bam:
    input:
        get_merge_bam_inputs,
    wildcard_constraints:
        sample="|".join(samples.keys())  # Only match actual sample names (not sample_run names)
    output:
        bam = os.path.join(result_path,"results","{sample}","mapped", "{sample}.filtered.bam"),
        bai = os.path.join(result_path,"results","{sample}","mapped", "{sample}.filtered.bam.bai"),
        stats = os.path.join(result_path, 'results', "{sample}", '{sample}.align.stats.tsv'),
    params:
        sample_name = "{sample}",
        mitochondria_name = config["mitochondria_name"],
        num_inputs = lambda w: len(get_runs_for_sample(w.sample)),
        input_bams = lambda w: " ".join([os.path.join(result_path, "results", sr, "mapped", f"{sr}.filtered.bam") 
                                          for sr in get_runs_for_sample(w.sample)]),
        # Only match if this is a valid sample name (check happens in input function)
        is_valid_sample = lambda w: w.sample in samples,
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 2)
    conda:
        "../envs/bowtie2.yaml",
    log:
        "logs/rules/merge_bam_{sample}.log"
    shell:
        """
        mkdir -p $(dirname {output.bam})
        
        if [ {params.num_inputs} -eq 1 ]; then
            # Single run, just copy
            cp {input[0]} {output.bam}
            if [ -f {input[0]}.bai ]; then
                cp {input[0]}.bai {output.bai}
            else
                samtools index {output.bam}
            fi
        else
            # Multiple runs, merge
            samtools merge -f -@ {threads} {output.bam} {params.input_bams} 2> {log}
            samtools index {output.bam} 2>> {log}
        fi
        
        # Calculate stats from merged BAM
        samtools idxstats {output.bam} | awk '{{ sum += $3 + $4; if($1 == "{params.mitochondria_name}") {{ mito_count = $3; }}}}END{{ print "mitochondrial_fraction\t"mito_count/sum }}' > {output.stats}
        """
        
# peak calling with MACS2 & samtools and annotation with HOMER
rule peak_calling:
    input:
        bam = os.path.join(result_path,"results","{sample}","mapped", "{sample}.filtered.bam"),
        bai = os.path.join(result_path,"results","{sample}","mapped", "{sample}.filtered.bam.bai"),
        homer_script = os.path.join(HOMER_path,"configureHomer.pl"),
        regulatory_regions = config["regulatory_regions"],
    output:
        peak_calls = os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.narrowPeak"),
        peak_annot = os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.narrowPeak.annotated.tsv"),
        peak_annot_log = os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.narrowPeak.annotated.tsv.log"),
        macs2_xls = os.path.join(result_path,"results","{sample}","peaks","{sample}_peaks.xls"),
        summits_bed = os.path.join(result_path,"results","{sample}","peaks","{sample}_summits.bed"),
        homer_knownResults = os.path.join(result_path,"results","{sample}","homer","knownResults.txt"),
        homer_log = os.path.join(result_path,"results","{sample}","homer","{sample}.homer.log"),
        macs2_log = os.path.join(result_path, 'results', "{sample}", 'peaks', '{sample}.macs2.log'),
        stats = os.path.join(result_path, 'results', "{sample}", '{sample}.peak.stats.tsv'),
    params:
        peaks_dir = os.path.join(result_path,"results","{sample}","peaks"),
        homer_dir = os.path.join(result_path,"results","{sample}","homer"),
        homer_bin = os.path.join(HOMER_path,"bin"),
        formating = lambda w: '--format BAMPE' if samples["{}".format(w.sample)]["read_type"] == "paired" else '--format BAM',
        genome_size = config["genome_size"],
        genome = config["genome"],
        keep_dup = config['macs2_keep_dup'],
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: 4*config.get("threads", 2)
    conda:
        "../envs/macs2_homer.yaml",
    log:
        "logs/rules/peak_calling_{sample}.log"
    shell:
        """
        export PATH="{params.homer_bin}:$PATH";
        
        macs2 callpeak -t {input.bam} {params.formating} \
            --nomodel --keep-dup {params.keep_dup} --extsize 147 -g {params.genome_size} \
            -n {wildcards.sample} \
            --outdir {params.peaks_dir} > "{output.macs2_log}" 2>&1;
        
        {params.homer_bin}/annotatePeaks.pl {output.peak_calls} {params.genome} > {output.peak_annot} 2> {output.peak_annot_log};
        
        {params.homer_bin}/findMotifsGenome.pl "{output.summits_bed}" {params.genome} {params.homer_dir} -size 200 -mask > "{output.homer_log}" 2>&1

        cat {output.peak_calls} | wc -l | awk '{{print "peaks\t" $1}}' >> "{output.stats}"
        
        TOTAL_READS=`samtools idxstats {input.bam} | awk '{{sum += $3}}END{{print sum}}'`;
        
        samtools view -c -L {output.peak_calls} {input.bam} | awk -v total=$TOTAL_READS '{{print "frip\t" $1/total}}' >> "{output.stats}";

        samtools view -c -L {input.regulatory_regions} {input.bam} | awk -v total=$TOTAL_READS '{{print "regulatory_fraction\t" $1/total}}' >> "{output.stats}";
        
        if [ ! -f {output.homer_knownResults} ]; then
            touch {output.homer_knownResults}
        fi
        """
        
rule aggregate_stats:
    input:
        align_stats = os.path.join(result_path, 'results', "{sample}", '{sample}.align.stats.tsv'),  # From merge_bam
        peak_stats = os.path.join(result_path, 'results', "{sample}", '{sample}.peak.stats.tsv'),
    output:
        os.path.join(result_path, 'results', "{sample}", '{sample}.stats.tsv'),
    resources:
        mem_mb=config.get("mem", "1000"),
    threads: config.get("threads", 2)
    log:
        "logs/rules/aggregate_stats_{sample}.log"
    shell:
        """
        cat {input.align_stats} {input.peak_stats} > {output}
        """
