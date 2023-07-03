#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// set up and create an output directory
outdir = file(params.outdir)
outdir.mkdir()

log.info """\

        ================================================================
        Rapid CNS2 - Nextflow P I P E L I N E
        ================================================================

        INPUTS
        ================================================================
        sample        : ${params.sample}
        input_bam     : ${params.bam}
        variant_vcf   : ${params.variant_vcf}
        bed_methyl    : ${params.bedmethyl}
        outdir        : ${params.outdir}
        patient       : ${params.patient}
        threads       : ${params.threads}
        ================================================================

        """
        .stripIndent()

process samtools_index {
    input:
        path(params.bam)

    output:
        path "*.bai", emit: indexed_bam

    script:
        """
        samtools index ${params.bam}
        """
}

process cnvpytor{
    input:
        val(sample)
        path(input_bam)
        val(threads)

    output:
        val true
    
    publishDir("${params.outdir}")

    script:
        """
        cnvpytor -root ${sample}_CNV.pytor -rd ${input_bam.toRealPath()} -j ${threads}
        cnvpytor -root ${sample}_CNV.pytor -his 1000 10000 100000 -j ${threads} # 4 mins
        cnvpytor -root ${sample}_CNV.pytor -partition 1000 10000 100000 -j ${threads} # SLOW
        cnvpytor -root ${sample}_CNV.pytor -call 1000 -j ${threads} > ${sample}.cnvpytor.calls.1000.tsv
        cnvpytor -root ${sample}_CNV.pytor -call 10000 -j ${threads} > ${sample}.cnvpytor.calls.10000.tsv
        cnvpytor -root ${sample}_CNV.pytor -call 100000 -j ${threads} > ${sample}.cnvpytor.calls.100000.tsv
        cnvpytor -root ${sample}_CNV.pytor -plot manhattan 100000 -chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 -o ${PWD}/${params.outdir}/${sample}_cnvpytor_100k.png
        """
}

process vcftools {
    input:
        path(variant_vcf)
        val(sample)

    output:
        path "*.vcf", emit: variants_out

    script:
        """
        vcftools --gzvcf ${variant_vcf} --remove-filtered-all --recode --out ${sample}
        """
}

process gzip {
    input:
        path(input_file)

    output:
        path "*.gz", emit: compressed_out

    script:
        """
        gzip -c ${input_file} > ${input_file}.gz
        """
}

process bedtools_intersect {
    input:
        path(input1)
        path(input2)
        val(output_file)
        val(ext)

    output:   
        path "*.vcf", emit: intersect_vcf
        
    script:
        """
        bedtools intersect -a ${input1} -b ${input2} > ${output_file}${ext}
        """
}

process bedtools_intersect2 {
    input:
        path(input1)
        path(input2)
        val(output_file)
        val(ext)

    output:   
        path "*.bed", emit: intersect_bed
        
    script:
        """
        bedtools intersect -a ${input1} -b ${input2} > ${output_file}.${ext}
        """
}

process convert2annovar{
    input:
        path(input)
        val(output_file)
        val(ext)

    output:
        path "*.avinput", emit: annovar_input

    script:
        """
        /annovar/convert2annovar.pl -format vcf4 ${input} -withfreq -includeinfo > ${output_file}${ext}
        """
}

process table_annovar {
    input:
        path(input)
        val(annovar_ver)
        val(output_file)
        val(ext)
    
    output:
        path "*_multianno.csv", emit: clair3_output
      
    script:
        """
        /annovar/table_annovar.pl ${input} /annovar/humandb/ -buildver ${annovar_ver} -out ${output_file}${ext} -protocol refGene,cytoBand,avsnp147,dbnsfp30a,1000g2015aug_eur,cosmic70 -operation gx,r,f,f,f,f -nastring . -csvout -polish -otherinfo
        """
}

process mgmt_pred {
    input:
        path(mgmt_pred)
        file(intersect_bed)
        path(probes)
        path(model)
        val(sample)
        val(params.outdir)
    
    publishDir("${params.outdir}")

    output:
        val true

    script:
        """
        Rscript ${mgmt_pred} --input ${intersect_bed} --probes ${probes} --model ${model} --out_dir ${PWD}/${params.outdir} --sample ${sample} 
        """
}

process meth_classification {
    input:
        path(meth_class)
        val(sample)
        val(params.outdir)
        path(bedmethyl)
        path(topprobes)
        path(trainingdata)
        path(arrayfile)
        val(threads)

    publishDir("${params.outdir}")

    output:
        val true

    script:
        """
        Rscript ${meth_class} --sample ${sample} --out_dir ${PWD}/${params.outdir} --in_file ${bedmethyl} --probes ${topprobes} --training_data ${trainingdata} --array_file ${arrayfile} --threads ${threads}
        """
}

process filter_report {
    input:
        path(filterreport)
        path(clair3_multianno)
        val(sample)
        val(params.outdir)
    
    publishDir("${params.outdir}")

    output:
        val true
    
    script:
        """
        Rscript ${filterreport} --input ${clair3_multianno} --output ${PWD}/${params.outdir}/${sample}_clair3_report.csv --sample ${sample}
        """        
}

process make_report {
    input:
        path(makereport)
        val(ready) // cnvpytor
        val(ready) // mgmt_pred
        val(ready) // meth_classification
        val(ready) // filter_report
        val(sample)
        path(outdir)
        path(report_UKHD)
    
    output:

    publishDir("${params.outdir}")

    script:
        """
        Rscript ${makereport} \
            --prefix ${sample} \
            --mutations ${PWD}/${params.outdir}/${sample}_clair3_report.csv \
            --cnv_plot ${PWD}/${params.outdir}/${sample}_cnvpytor_100k.global.0000.png \
            --rf_details ${PWD}/${params.outdir}/${sample}_rf_details.tsv \
            --votes ${PWD}/${params.outdir}/${sample}_votes.tsv \
            --output_dir ${PWD}/${params.outdir} \
            --patient ${sample} \
            --coverage ${PWD}/${params.outdir}/${sample}.mosdepth.summary.txt \
            --sample ${sample} \
            --mgmt ${PWD}/${params.outdir}/${sample}_mgmt_status.csv \
            --report_UKHD ${report_UKHD}
        """
}

///////////////////////////
// MAIN WORKFLOW SECTION //
///////////////////////////

workflow {
/////////////////////
    // input(s) 
    Channel.fromPath(params.outdir)
    .set {outdir}

    Channel.fromPath("${projectDir}/bin/NPHD_panel_hg38_clean.bed", checkIfExists: true)
    .set {targets}

    Channel.fromPath(params.bam, checkIfExists: true)
    .set {input_bam}

    Channel.from(params.sample)
    .set {sample}

    Channel.fromPath(params.mosdepth_data, checkIfExists: true)
    .set {mosdepth_data}

    Channel.fromPath(params.variant_vcf, checkIfExists: true)
    .set {variant_vcf}

    Channel.fromPath(params.bedmethyl, checkIfExists: true)
    .set {bedmethyl}

    Channel.fromPath("${projectDir}/bin/mgmt_hg38.bed", checkIfExists: true)
    .set {mgmt_bed}

    Channel.fromPath("${projectDir}/bin/mgmt_pred_v0.1.R", checkIfExists: true)
    .set {mgmt_pred}

    Channel.fromPath("${projectDir}/bin/mgmt_probes.Rdata", checkIfExists: true)
    .set {probes}

    Channel.fromPath("${projectDir}/bin/mgmt_137sites_mean_model.Rdata", checkIfExists: true)
    .set {model}
    
    Channel.fromPath("${projectDir}/bin/methylation_classification_nanodx_v0.1.R", checkIfExists: true)
    .set {meth_class}

    Channel.fromPath("${projectDir}/bin/top_probes_hm450.Rdata", checkIfExists: true)
    .set {topprobes}
    
    Channel.fromPath("${projectDir}/bin/capper_top_100k_betas_binarised.Rdata", checkIfExists: true)
    .set {trainingdata}

    Channel.fromPath("${projectDir}/bin/HM450.hg38.manifest.gencode.v22.Rdata", checkIfExists: true)
    .set {arrayfile}

    Channel.from(params.threads)
    .set {threads}
    
    Channel.fromPath("${projectDir}/bin/filter_report_v0.1.R", checkIfExists: true)
    .set {filterreport}

    Channel.fromPath("${projectDir}/bin/make_report_v0.1.R", checkIfExists: true)
    .set {makereport}

    Channel.from(params.patient)
    .set {patient}

    Channel.fromPath("${projectDir}/bin/Rapid_CNS2_report_UKHD_v0.1.Rmd", checkIfExists: true)
    .set {report_UKHD}

///////////////////// - run the workflow

    // index the input bam
    // NOTE - don't be tempted to remove this (if bam already indexed)
    // removing this indexing step will break CNVpytor    
    samtools_index_ch = samtools_index(input_bam)

    // run cnvPYTOR
    cnvpytor(sample, input_bam, threads)

    // run annovar
    variants_ch = vcftools(variant_vcf, sample)

    // compress output from annovar
    compressed_variants_ch = gzip(variants_ch.variants_out)

    // run bedtools intersect on output of vcftools
    vcf_intersect_ch = bedtools_intersect(compressed_variants_ch.compressed_out, targets, sample, '_clair3_panel.vcf')
    
    // convert vcf to annovar input
    converted_annovar_ch = convert2annovar(vcf_intersect_ch.intersect_vcf, sample, '_clair3_panel.avinput')

    // run table_annovar on the converted annovar input
    clair3_annovar_ch = table_annovar(converted_annovar_ch.annovar_input, 'hg38', sample, '_clair3_panel')

    // run bedtools intersect
    intersect_bed_ch = bedtools_intersect2(bedmethyl, mgmt_bed, 'mgmt_5mC.hg38', 'bed')

    // run the mgmt_pred script
    mgmt_pred(mgmt_pred, intersect_bed_ch.intersect_bed, probes, model, sample, params.outdir)

    // run the meth classification script
    meth_classification(meth_class, sample, params.outdir, bedmethyl, topprobes, trainingdata, arrayfile, threads)

    // collect report data and generate report
    filter_report(filterreport, clair3_annovar_ch.clair3_output, sample, params.outdir)  

    // collect data and generate final report
    make_report(makereport, cnvpytor.out, mgmt_pred.out, meth_classification.out, filter_report.out, sample, outdir, report_UKHD)
}