ROOT=/net/topmed11/working/porchard/variant-sensitive-motif-scanning
DATA=$(ROOT)/data
WORK=$(ROOT)/work
BIN=$(ROOT)/bin
PIPELINES=$(ROOT)/pipelines
ANALYSIS=$(WORK)/$@

define NL


endef

motifs:
	cd $(DATA) && tar -xvzf motifs.tar.gz

vcf:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --vcf_glob '/net/topmed10/working/porchard/rnaseq/work/scan-variant-vcf-files/freeze-2b/results/vcfs-by-chrom/Whole_blood*.vcf.gz' $(ROOT)/make-vcf.nf &

fasta:
	mkdir -p $(DATA)/fasta && ln -s /net/topmed10/working/porchard/rnaseq/data/fasta/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta $(DATA)/fasta/

scan:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume -qs 1000 --results $(ANALYSIS)/results --meme_glob '$(DATA)/motifs/*.meme' --fasta $(DATA)/fasta/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta --vcf $(DATA)/vcf/*.vcf.gz $(PIPELINES)/snp-sensitive-motif-scanning/main.nf &

fimo-ref-and-alt-scores:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume -qs 1000 --results $(ANALYSIS)/results --fimo_ref_glob '$(WORK)/scan/results/scan/ref/*.txt' --fimo_alt_glob '$(WORK)/scan/results/scan/alt/*.txt' --vcf $(WORK)/vcf/results/bcfs-merged/merged.vcf.gz $(ROOT)/get-fimo-ref-and-alt-scores.nf &
