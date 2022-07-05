include { bam_reset } from '../modules/bam_reset.nf'


workflow reset_bam_alignment {

	//remove alignment from bam - this process proceeds directly after the end of 1.2

	take:
		input_manifest //from step 1.2, path/file object

	main:
	
		//declare input_ch to 1.3
		input_ch = Channel.fromPath( input_manifest, checkIfExists : true )
	       	       .splitCsv( header : true )
	               .map { 
				row -> sample_tag: row.sample_tag
				       bam_ch: row.bam_fl
			                         
			 }

		sample_tag = input_ch.map { it -> it[0] }
		bam_ch = input_ch.map { it -> it[1] }

	//reset bam alignment

		bam_reset(sample_tag, bam_ch)			   

}

workflow {


        manifest_ch = Channel.fromPath("./test_manifest.csv")

        reset_bam_alignment(manifest_ch)


}

