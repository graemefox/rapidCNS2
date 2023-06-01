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
        variant_vcf   : ${params.variant_vcf}
        bed_methyl    : ${params.bedmethyl}
        outdir        : ${params.outdir}
        patient       : ${params.patient}
        cnvpytor_plot : ${params.cnvpytor_plot}
        mosdepth_plot : ${params.mosdepth_plot}
        threads       : ${params.threads}
        ================================================================

        """
        .stripIndent()

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
    
    publishDir("${params.outdir}")

    output:
        path "*_multianno.csv", emit: clair3_output
      
    script:
        """
        /annovar/table_annovar.pl ${input} /annovar/humandb/ -buildver ${annovar_ver} -out ${output_file}${ext} -protocol refGene,cytoBand,avsnp147,dbnsfp30a,1000g2015aug_eur,cosmic70 -operation gx,r,f,f,f,f -nastring . -csvout -polish -otherinfo
        """
}

process R_mgmt_pred {
    input:
        path(mgmt_pred)
        file(intersect_bed)
        path(probes)
        path(model)
        val(sample)
        path(params.outdir)
    
    publishDir("${params.outdir}")

    output:
        val '*mgmt_status.csv', emit: mgmt_status

    script:
        """
        Rscript ${mgmt_pred} --input ${intersect_bed} --probes ${probes} --model ${model} --out_dir ${params.outdir} --sample ${sample} 
        """
}

process R_meth_classification {
    input:
        path(meth_class)
        val(sample)
        path(outdir)
        path(bedmethyl)
        path(topprobes)
        path(trainingdata)
        path(arrayfile)
        val(threads)

    publishDir("${params.outdir}")

    output:
        val "*_rf_details.tsv", emit: rf_details
        val "*_votes.tsv", emit: votes

    script:
        """
        Rscript ${meth_class} --sample ${sample} --out_dir ${params.outdir} --in_file ${bedmethyl} --probes ${topprobes} --training_data ${trainingdata} --array_file ${arrayfile} --threads ${threads}
        """
}

process filter_report {
    input:
        path(filterreport)
        path(clair3_multianno)
        val(sample)
    
    publishDir("${params.outdir}")

    output:
        val "*_clair3_report.csv", emit: clair3_report
    
    script:
        """
        Rscript ${filterreport} --input ${clair3_multianno} --output ${sample}_clair3_report.csv --sample ${sample}
        """        
}

process make_report {
    input:
        path(makereport)
        path(report_UKHD)
        val(sample)
        path(params.outdir)
        val(mgmt_status)
        val(clair3_report)
        val(rfdetails)
        val(votes)
        path(cnvpytor_plot)
        path(mosdepth_plot)
        val(patient)

    output:

    publishDir("${params.outdir}")

    script:
        """
        Rscript ${makereport} --rf_details=${params.outdir}/${sample}_rf_details.tsv --votes=${params.outdir}/${sample}_votes.tsv --mutations=${params.outdir}/${sample}_clair3_report.csv --coverage=${mosdepth_plot} --output_dir=${params.outdir}/${sample}_report/ --cnv_plot=${cnvpytor_plot} --prefix=${sample} --mgmt=${params.outdir}/${sample}_mgmt_status.csv --patient=${patient} --sample=${sample} --report_UKHD=${report_UKHD}
        """
        //
        //Rscript ${makereport} --rf_details=${rfdetails} --votes=${votes} --mutations=${clair3_report} --coverage=${mosdepth_plot} --output_dir=${params.outdir}/${sample}_report/ --cnv_plot=${cnvpytor_plot} --prefix=${sample} --mgmt=${mgmt_status} --patient=${patient} --sample=${sample} --report_UKHD=${report_UKHD}
}

///////////////////////////
// MAIN WORKFLOW SECTION //
///////////////////////////

workflow {

/////////////////////  -  collect variables
    // input(s) for annovar process
    Channel.fromPath(params.variant_vcf, checkIfExists: true)
    .set {variant_vcf}

    Channel.from(params.sample)
    .set {sample}

    Channel.fromPath(params.outdir)
    .set {outdir}

    // TODO  -  check this 'CLEAN' panel - is this what we want????
    Channel.fromPath("${projectDir}/bin/NPHD_panel_hg38_clean.bed", checkIfExists: true)
    .set {targets}

/////////////////////
    // input(s) for bedtools_intersect process
    Channel.fromPath(params.bedmethyl, checkIfExists: true)
    .set {bedmethyl}

    Channel.fromPath("${projectDir}/bin/mgmt_hg38.bed", checkIfExists: true)
    .set {mgmt_bed}

/////////////////////
    // input(s) for R_mgmt_pred process
    Channel.fromPath("${projectDir}/bin/mgmt_pred_v0.1.R", checkIfExists: true)
    .set {mgmt_pred}

    Channel.fromPath("${projectDir}/bin/mgmt_probes.Rdata", checkIfExists: true)
    .set {probes}

    Channel.fromPath("${projectDir}/bin/mgmt_137sites_mean_model.Rdata", checkIfExists: true)
    .set {model}
    
/////////////////////
    // input(s) for R_meth_classification process
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
    
/////////////////////
    // Inputs for generating reports
    Channel.fromPath("${projectDir}/bin/filter_report_v0.1.R", checkIfExists: true)
    .set {filterreport}

    Channel.fromPath("${projectDir}/bin/make_report_v0.1.R", checkIfExists: true)
    .set {makereport}

    Channel.from(params.patient, checkIfExists: true)
    .set {patient}

    Channel.fromPath("${projectDir}/bin/Rapid_CNS2_report_UKHD_v0.1.Rmd", checkIfExists: true)
    .set {report_UKHD}

    Channel.fromPath(params.cnvpytor_plot, checkIfExists: true)
    .set {cnvpytor_plot}

    Channel.fromPath(params.mosdepth_plot, checkIfExists: true)
    .set {mosdepth_plot}

///////////////////// - run the workflow
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
    mgmt_pred_ch = R_mgmt_pred(mgmt_pred, intersect_bed_ch.intersect_bed, probes, model, sample, params.outdir)

    // run the meth classification script
    meth_class_ch = R_meth_classification(meth_class, sample, params.outdir, bedmethyl, topprobes, trainingdata, arrayfile, threads)

    // collect report data and generate report
    filter_report_ch = filter_report(filterreport, clair3_annovar_ch.clair3_output, sample)  

    // NOTE - the channel inputs aren't used (report_UKHD not happy with inputs being passed as channels)
    // but they delay this being run until the inputs have been generated
    //make_report(makereport, report_UKHD, sample, params.outdir, mgmt_pred_ch.mgmt_status, filter_report_ch.clair3_report, meth_class_ch.rf_details, meth_class_ch.votes, cnvpytor_plot, mosdepth_plot)  
    make_report(makereport, report_UKHD, sample, params.outdir, mgmt_pred_ch.mgmt_status, filter_report_ch.clair3_report, meth_class_ch.rf_details, meth_class_ch.votes, cnvpytor_plot, mosdepth_plot, patient)  
}