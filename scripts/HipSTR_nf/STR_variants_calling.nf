
/*
~ ~ ~ > * USER INPUT PARAMETERS 
*/
date = new Date().format( 'yyyyMMdd' )

params.bamfiles  	= null
params.ref_fa 		= "bin/c_elegans.PRJNA13758.WS276.genomic.fa"
params.ref_str 		= "build_ref/ref_ce.hipstr_reference.bed"
params.email = ""

/*
~ ~ ~ > * OUTPUT DIRECTORY 
*/

params.out = "STR_Results-${date}"




log.info ""
log.info "------------------------------------------"
log.info "        STR variant calling pipeline "
log.info "------------------------------------------"
log.info ""



log.info ""
log.info "bam files                          = ${params.bamfiles}"
log.info "Genome fasta                       = ${params.ref_fa}"
log.info "Reference STR                      = ${params.ref_str}"
log.info ""










/*
~ ~ ~ > * INITIATE bam channel
*/

File bams = new File("${params.bamfiles}")
bam_handle = bams.getAbsolutePath()



/*
~ ~ ~ > * INITIATE ref fasta channel
*/

File reffa = new File("${params.ref_fa}")
reffa_handle = reffa.getAbsolutePath()


/*
~ ~ ~ > * INITIATE ref STR channel
*/

File refstr = new File("${params.ref_str}")
refstr_handle = refstr.getAbsolutePath()




CONTIG_LIST = ["I", "II", "III", "IV", "V", "X"]

contig = Channel.from(CONTIG_LIST)




process split_bed {


  cpus 1
  memory '2 GB'




  output: 

   file("*.00*.bed") into bed_files

  """



  cat ${refstr_handle} | grep -v "MtDNA"   > refstr_noMt.bed

  for chr in `cut -f 1 refstr_noMt.bed | sort | uniq`; do
                grep -w \$chr refstr_noMt.bed > refstr_noMt_\$chr.bed

  bedtools split -i refstr_noMt_\$chr.bed -n 20 -p \$chr


  done


 
  """
}







bed_files
.flatten()
.map { file -> tuple(file.baseName.replaceAll(/.bed/,""), file, file.getSimpleName() ) }
.into{chrbed; chrbed_print}

 

process str_call {


	cpus 4
	memory '32 GB'

	input:

	set val(range),  file(refstr_part), val(CHROM) from chrbed

	output: 

	file "str_${CHROM}_${range}.vcf.gz" into str_6chr


	"""

  

  HipSTR --bam-files ${bam_handle} \\
       --fasta ${reffa_handle} \\
       --regions ${refstr_part} \\
       --output-filters \\
       --chrom ${CHROM} \\
       --str-vcf str_${CHROM}_${range}.vcf.gz \\
       --log log_${CHROM}_${range}.txt \\
       --stutter-out stutter_models_${CHROM}_${range}.txt \\
       --viz-out viz_${CHROM}_${range}.viz.gz



	"""
}





process MtDNA_str_call {

  cpus 8
  memory '64 GB'

  output: 

  file "str_MtDNA.vcf.gz" into str_Mt

  """



cat ${refstr_handle} | grep "MtDNA" > hipstr_ref_MtDNA.bed


  HipSTR --bam-files ${bam_handle} \\
       --fasta ${reffa_handle} \\
       --regions hipstr_ref_MtDNA.bed \\
       --output-filters \\
       --str-vcf str_MtDNA.vcf.gz \\
       --log log_MtDNA.txt \\
       --stutter-out stutter_models_MtDNA.txt \\
       --viz-out viz_MtDNA.viz.gz
  """
}




process merge_filter {

  publishDir "${params.out}", mode: 'copy', pattern: "STR_all_*.vcf.*"


  cpus 4
  memory '32 GB'


  input:

  file('*') from str_6chr.mix(str_Mt).collect()

  output: 

  file "STR_all_*.vcf.gz" 

  """



  bcftools concat str_*.vcf.gz | bcftools sort -Oz -o STR_all_raw.vcf.gz

  python $PWD/bin/HipSTR_filter_vcf.py --vcf STR_all_raw.vcf.gz --min-call-qual 0.9 --max-call-flank-indel 0.15 --max-call-stutter 0.15 --min-call-allele-bias -2 --min-call-strand-bias -2 > STR_filtered.vcf.gz


  bcftools view STR_filtered.vcf.gz -Oz -o STR_all_filtered.vcf.gz


  bcftools view STR_all_filtered.vcf.gz | bcftools filter -i "F_MISSING<0.1" -Oz -o STR_all_filtered_Fmiss01.vcf.gz
 
 

  """
}











workflow.onComplete {

    summary = """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """

    println summary


    // mail summary
    if (params.email) {
        ['mail', '-s', 'str-nf', params.email].execute() << summary
    }//


}
