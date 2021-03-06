The Leptospira variant calling pipeline is a bash pipeline that use state of the art tools to obtain small variants (SNVs and short indels) and Copy Number Variants (CNVs) from Leptospira NGS DNA sequencing data.
The input format is BAM.
The output formats are VCF for small variants and GFF for CNVs.

## Dependencies

**Required**
* Oracle JDK 1.8
* GATK 3.5
* Picard 2.2.2
* SnpEff 4.2
* Samtools 1.3.1
* Bcftools 1.3.1

**Optional**
* FastQC v0.11.5
* CNVnator v0.3 (required for CNV calling)

## Installation

Intall the dependencies.
Clone the repository. Before running the pipeline, update the paths to dependencies in `config/config.sh` and source it:

	source config/config.sh

Consider adding this line to `~/.bashrc`.

## Running the pipeline

### SNVs and short indel variant calling

To call for SNPs and short indels either run the `scvc.sh` script:

	scvc.sh <INPUT_BAM> <OUTPUT_DIR> <REFERENCE>

or alternatively run each step separately:

1) Prepare the reference (creates all required indices, input file must have .fasta extension): `reference/prepare_reference.sh <INPUT_FASTA>`

2) Prepare BAM files (clean, fix mate info, remove duplicates and then realigns it around indels): `prepare_bam.sh <INPUT_BAM> <OUTPUT_FOLDER> <REFERENCE>`

3) Run the variant callers on the realigned BAM (realignments\INPUT.bam.realigned.bam):
```
samtools_pileup.sh <INPUT_BAM> <ST_OUTPUT_VCF> <REFERENCE>
unified_genotyper.sh <INPUT_BAM> <HC_OUTPUT_VCF> <REFERENCE>
haplotype_caller.sh <INPUT_BAM> <UG_OUTPUT_VCF> <REFERENCE>
```

4) Filter the VCF files using recommended hard thresholds:
```
variant_filtering_samtools.sh <ST_INPUT_VCF> <ST_OUTPUT_VCF> <REFERENCE>
variant_filtering_gatk.sh <HC_INPUT_VCF> <HC_OUTPUT_VCF> <REFERENCE>
variant_filtering_gatk.sh <UG_INPUT_VCF> <UG_OUTPUT_VCF> <REFERENCE>
```

5) Combine variant calls from HaplotypeCaller, UnifiedGenotyper and samtools in the intersection and union sets: `combine_variants.sh <HC_INPUT_VCF> <UG_INPUT_VCF> <ST_INPUT_VCF> <UNION_OUTPUT_VCF> <INTERSECTION_OUTPUT_VCF> <REFERENCE>`

6) Run SNPEff for variant annotation: `annotation.sh <INPUT_VCF> <OUTPUT_VCF> <SNPEFF_REFERENCE> <REFERENCE>`, where alternatives for `<SNPEFF_REFERENCE>` are `Leptospira_borgpetersenii_serovar_Hardjo_bovis_L550_uid58507` or `Leptospira_borgpetersenii_serovar_Hardjo_bovis_JB197_uid58509`.

For example, to run the entire pipeline for a `LBH-A_JB197.bam` file resulting from aligning the reads of a strain LBH-A against reference genome `Lb.Hardjo.JB197.fasta`:

	scvc.sh LBH-A_JB197.bam variants Lb.Hardjo.JB197.fasta

or to run each step separately:

	prepare_bam.sh LBH-A_JB197.bam realignments/ Lb.Hardjo.JB197.fasta > LBH-A_JB197_prepare_bam.log
	variant_calling/samtools_pileup.sh realignments/LBH-A_JB197.bam.realigned.bam variants/samtools/LBH-A_JB197.raw.vcf Lb.Hardjo.JB197.fasta > variants/samtools/LBH-A_JB197_samtools.log
	variant_calling/unified_genotyper.sh realignments/LBH-A_JB197.bam.realigned.bam variants/unified_genotyper/LBH-A_JB197.raw.vcf Lb.Hardjo.JB197.fasta > variants/unified_genotyper/LBH-A_JB197_ug.log
	variant_calling/haplotype_caller.sh realignments/LBH-A_JB197.bam.realigned.bam variants/haplotype_caller/LBH-A_JB197.raw.vcf Lb.Hardjo.JB197.fasta > variants/haplotype_caller/LBH-A_JB197_hc.log
	variant_filtering/variant_filtering_samtools.sh variants/samtools/LBH-A_JB197.raw.vcf variants/samtools/LBH-A_JB197.filtered.vcf Lb.Hardjo.JB197.fasta
	variant_filtering/variant_filtering_gatk.sh variants/haplotype_caller/LBH-A_JB197.raw.vcf variants/haplotype_caller/LBH-A_JB197.filtered.vcf Lb.Hardjo.JB197.fasta
	variant_filtering/variant_filtering_gatk.sh variants/unified_genotyper/LBH-A_JB197.raw.vcf variants/unified_genotyper/LBH-A_JB197.filtered.vcf Lb.Hardjo.JB197.fasta
	combine_variants/combine_variants.sh variants/haplotype_caller/LBH-A_JB197.filtered.vcf variants/unified_genotyper/LBH-A_JB197.filtered.vcf variants/samtools/LBH-A_JB197.filtered.vcf variants/LBH-A_JB197.union.vcf variants/LBH-A_JB197.intersection.vcf Lb.Hardjo.JB197.fasta
	annotation/annotation.sh variants/LBH-A_JB197.union.vcf variants/LBH-A_JB197.union.annotated.vcf Leptospira_borgpetersenii_serovar_Hardjo_bovis_JB197_uid58509 Lb.Hardjo.JB197.fasta

### CNV variant calling

To call for CNVs using CNVnator:

	cnvnator.py [-h] [--window_size <n>] <input_bam> <input_reference> <output_folder>

This runs the CNVnator CNV calling pipeline. Outputs CNVs in CNVnator native format and GFF format.

Positional arguments:

* input_bam:             Input BAM alignments file
* input_reference:       Input FASTA reference file
* output_folder:         Output folder

Optional arguments:

* -h, --help:           show this help message and exit
* --window_size <n>     The window size should be determined by the average
                        read depth and the read length. The recommended values
                        are as follows ~100-bp for 20-30x coverage, ~500-bp
                        for 4-6x coverage, and ~30-bp bins for 100x coverage.
                        This value should not be lower than the read length
                        though.

Example:

`python src/hcvc/cnv_calling/cnvnator.py ~/data/BK-30_L550.bam.realigned.bam ~/data/Lb.Hardjo.L550.fasta ~/data/output_folder --window_size 300`


## Methodology

### Small variant filtering

There are different filtering criteria for different variant callers and for different types of variants. Variant filtering is not implemented for CNVs.

SNVs from GATK UnifiedGenotyper and Haplotype Caller are filtered under any of the following conditions:
* `QD < 2.0`
* `SOR > 6.0`
* `QUAL < 5`
* `DP < 3`

Indels from GATK UnifiedGenotyper and Haplotype Caller are filtered under any of the following conditions:
* `QD < 2.0`
* `SOR > 10.0`
* `QUAL < 5`
* `DP < 3`

SNVs from Samtools are filtered under any of the following conditions:
* `QUAL < 17`
* `DP < 3`

Indels from Samtools are filtered under any of the following conditions:
* `QUAL < 50`
* `DP < 3`

`QUAL` is the Phred quality score for the variant calling.
`DP` represents the depth of coverage.
For other technical annotations see GATK documentation [https://software.broadinstitute.org/gatk/documentation/tooldocs/3.5-0/](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.5-0/)

### Merging results from multiple variant callers

The pipeline outputs both the union and the intersection of the three variant callers. Intersection is usually a very restrictive alternative, the union has a higher sensitivity with a reasonable specificity.