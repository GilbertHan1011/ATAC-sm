import yaml

# Test the aligner switch logic
config_path = "config/config.yaml"

with open(config_path, 'r') as f:
    config = yaml.safe_load(f)

aligner = config.get("aligner", "bowtie2")
bwa_mem2_path = config.get("bwa_mem2_path")

print(f"Current aligner setting: {aligner}")
print(f"BWA-MEM2 index path: {bwa_mem2_path}")
print()

# Test bowtie2 rule
print("Testing rule 'align' (Bowtie2):")
if aligner == "bowtie2":
    print("  ✓ Rule 'align' should execute")
    print(f"  ✓ bowtie2_index input: {config['bowtie2_index']}")
else:
    print(f"  ✗ Rule 'align' will exit with error (aligner={aligner})")
    print("  ✗ bowtie2_index input: __DISABLED__/bowtie2_index (non-existent)")

print()

# Test bwa-mem2 rule
print("Testing rule 'align_bwa_mem' (BWA-MEM2):")
if aligner == "bwa-mem2":
    print("  ✓ Rule 'align_bwa_mem' should execute")
    if bwa_mem2_path and bwa_mem2_path != "null" and bwa_mem2_path != "":
        print(f"  ✓ Will use bwa-mem2 binary with index: {bwa_mem2_path}")
        print("  ✓ Index files: Not required (pre-indexed)")
    else:
        print(f"  ✓ Will use bwa binary with index: {config['genome_fasta']}")
        print("  ✓ Index files: Required from bwa_index rule")
else:
    print(f"  ✗ Rule 'align_bwa_mem' will exit with error (aligner={aligner})")

print()
print("Summary:")
print(f"  Active aligner: {aligner}")
if aligner == "bowtie2":
    print(f"  Active rule: align")
    print(f"  Binary: bowtie2")
elif aligner == "bwa-mem2":
    print(f"  Active rule: align_bwa_mem")
    if bwa_mem2_path and bwa_mem2_path != "null" and bwa_mem2_path != "":
        print(f"  Binary: bwa-mem2")
    else:
        print(f"  Binary: bwa")
