
rule fastqc_1:
    input:
        lambda w: get_units_fastqs(w)[0],
    output:
        html=os.path.join(result_path,"report","{sample}_fastqc_1.html"),
        zip=os.path.join(result_path,"report","{sample}_fastqc_1.zip"), # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "logs/rules/fastqc/{sample}.log",
    threads: 2
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.6.1/bio/fastqc"

rule fastqc_2:
    input:
        lambda w: get_units_fastqs(w)[1],
    output:
        html=os.path.join(result_path,"report","{sample}_fastqc_2.html"),
        zip=os.path.join(result_path,"report","{sample}_fastqc_2.zip"), # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "logs/rules/fastqc/{sample}.log",
    threads: 2
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.6.1/bio/fastqc"

rule trim_galore_pe:
    input:
        get_units_fastqs,
    output:
        fasta_fwd=os.path.join(result_path,"trimmed","{sample}_1.fq.gz"),
        report_fwd=os.path.join(result_path,"trimmed","{sample}_1._trimming_report.txt"),
        fasta_rev=os.path.join(result_path,"trimmed","{sample}_2.fq.gz"),
        report_rev=os.path.join(result_path,"trimmed","{sample}_2._trimming_report.txt"),
    threads: 2
    log:
        "logs/trim_galore/{sample}.log"
    wrapper:
        "v7.6.0/bio/trim_galore/pe"
