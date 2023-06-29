### Nextflow Rapid CNS2 pipeline

### Requirements:
```
docker
nextflow
```

### Pull docker image:
```
docker pull graefox/rapid_cns2:latest
```

### Pull latest workflow
```
nextflow pull graemefox/rapidCNS2
```


### Example commands:
```
## define sample name, ID and output directory:
SAMPLE=sample_01
PATIENT=JohnDoe
OUTPUT_DIR=/path/to/directory/${SAMPLE}_output

## run the pipeline
nextflow run graemefox/rapidCNS2 \
-with-docker graefox/rapid_cns2:latest \
-with-report ${OUTPUT_DIR}/${SAMPLE}_nextflow_report.html \
--bedmethyl /path/to/sample_01.methyl.cpg.bed.gz \
--variant_vcf /path/to/sample_01.wf_snp.vcf.gz \
--cnvpytor_plot /path/to/sample_01._cnvpytor_100k.global.0000.png \
--mosdepth_plot /path/to/sample_01.mosdepth.summary.txt \
--outdir ${OUTPUT_DIR} \
--sample ${SAMPLE} \
--patient ${PATIENT}
```
