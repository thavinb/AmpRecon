#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

include { bam_reset } from '../modules/bam_reset.nf'
include { bam_to_fastq } from '../modules/bam_to_fastq.nf'
include { align_bam } from '../modules/align_bam.nf'
workflow reset_bam_alignment {
  //remove alignment from bam - this process proceeds directly after the end of 1.2x

  take:
    sample_tag
    bam_file
    run_id
    reference_fasta_new
    reference_idx_fls
  main:
    // Unmap the bam files (ubam)
    bam_reset(sample_tag, bam_file)
    /*
    bam_to_fastq_In_ch = bam_reset.out.multiMap{
                            it -> (it.sample_tag), path(it.prealigned_bam)
                          }
    bam_to_fastq_In_ch.view()
    */

    // convert ubams to fastqs
    bam_to_fastq(bam_reset.out.sample_tag,
                 bam_reset.out.reset_bam)

    // do new alignment
    align_bam(
      bam_to_fastq.out.sample_tag,
      bam_to_fastq.out.fastq,
      reference_fasta_new,
      reference_idx_fls,
      run_id
      )
  //emit:

}
/*
// -------------------------- DOCUMENTATION -----------------------------------
[1] Collated BAM files are reset to their prealigned state by removing all
     @SQ from header, all reads marked as unmapped, dropping non-primary
     alignments and sorting order to set to unknown.
[2] BAM files are converted to FASTQ format, before being aligned to a reference genome.

*/

/*
// Simons branch code RIP

#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

include { bam_reset } from '../modules/bam_reset.nf'
include { align_bam } from '../modules/align_bam.nf'
include {indexReferenceBWA} from "../modules/bwa_index.nf"
include { bam_to_fastq } from '../modules/bam_to_fastq.nf'

workflow reset_bam_alignment {

	//remove alignment from bam - this process proceeds directly after the end of 1.2

	take:
		input_manifest //from step 1.2, path/file object what about irods

	main:
		input_ch = input_manifest.splitCsv(header : true)
					 .multiMap {
							row  -> bam_fl:row.bam_fl
							        sample_tag:row.sample_tag
						}
		 bam_reset(input_ch.sample_tag, input_ch.bam_fl)

}


workflow read_alignment {

	take:
		reference_genome
		input_fastq_ch //this can be output from bam_to_fastq


	main:
		sample_tag = input_fastq_ch.map { it -> it[0] }
		input_fastq = input_fastq_ch.map { it -> it[1] }
		reference_idx = indexReferenceBWA(reference_genome) //could add this as input to the workflow
		align_bam(sample_tag, input_fastq, reference_genome, reference_idx)

}
*/
