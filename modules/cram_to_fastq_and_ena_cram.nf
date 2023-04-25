process cram_to_fastq_and_ena_cram {
    /*
    * Converts a CRAM file into a FASTQ file and an ENA submission ready CRAM file
    */

    input:
        tuple val(file_id), path(sample_cram), val(reference_fasta)

    output:
        tuple val(file_id), path("${output_fastq_file}"), val(reference_fasta), emit: fastq
        path("${output_cram_file}"), emit: cram

    script:
        base_name=sample_cram.baseName
        output_fastq_file="${base_name}.fastq"
        output_cram_file="${base_name}.for_ena.cram"

        """
        # 1. Collate read pairs
        bamcollate2 collate=1 \
            level=9 \
            inputformat="cram" \
            filename=${sample_cram} \
            O="${base_name}.collated.bam"

        # 2. Removal of pre-existing alignments
        bamreset \
            resetaux=0 \
            level=9 \
            verbose=0 \
            < "${base_name}.collated.bam" \
            > "${base_name}.reset.bam"
        
        rm "${base_name}.collated.bam"

        # 3. Discarding of adapter sequences
        bamadapterclip \
            verbose=1 \
            level=0 \
            < "${base_name}.reset.bam" \
            > "${base_name}.adapter_clipped.bam"

        rm "${base_name}.reset.bam"

        # 4. Conversion to FASTQ format for output
        bamtofastq \
            < "${base_name}.adapter_clipped.bam" \
            > "${output_fastq_file}"

        # 5. Collate read pairs and attach post ranks to reads
        bamcollate2 collate=3 \
            level=9 \
            inputformat="bam" \
            filename="${base_name}.adapter_clipped.bam" \
            O="${base_name}.collated_and_ranked.bam"

        # 6. Conversion to FASTQ format
        bamtofastq \
            < "${base_name}.collated_and_ranked.bam" \
            > "${base_name}.intermediate.fastq"

        rm "${base_name}.collated_and_ranked.bam"

        # 7. Alignment to reference genome
        bwa mem \
            -p \
            -Y \
            -K 100000000 \
            -t 1 \
            "${reference_fasta}" \
            "${base_name}.intermediate.fastq" \
            > "${base_name}.realigned.sam"

        rm "${base_name}.intermediate.fastq"

        # 8. Conversion of SAM to BAM
        scramble -0 -I sam -O bam "${base_name}.realigned.sam" "${base_name}.realigned.bam"

        rm "${base_name}.realigned.sam"

        # 9. Merge headers
        merge_headers.py "${base_name}.adapter_clipped.bam" \
            "${base_name}.realigned.bam" \
            ${reference_fasta} \
            > "${base_name}.merged_header.txt"

        # 10. Reheader BAM.
        samtools reheader \
            "${base_name}.merged_header.txt" \
            "${base_name}.realigned.bam" \
            > "${base_name}.reheadered.bam"

        rm "${base_name}.merged_header.txt"
        rm "${base_name}.realigned.bam"

        # 11. Create single ranks per read 
        bam12split verbose=0 level=0 \
            < "${base_name}.reheadered.bam" \
            > "${base_name}.split.bam"

        rm "${base_name}.reheadered.bam"

        # 12. Merge BAM files
        bam12auxmerge \
            level=0 \
            rankstrip=1 \
            ranksplit=0 \
            zztoname=0 \
            clipreinsert=1 \
            "${base_name}.adapter_clipped.bam" \
            < "${base_name}.split.bam" \
            > "${base_name}.auxmerge.bam"

        rm "${base_name}.adapter_clipped.bam"
        rm "${base_name}.split.bam"

        # 13. Bambi select
        bambi select \
            --compression-level=0 \
            --input "${base_name}.auxmerge.bam" \
            --output "${base_name}.selected.bam" \
            -m "${base_name}.metrics"

        rm "${base_name}.auxmerge.bam"

        # 14. Coordinate sorting
        bamsort \
            threads=1 \
            SO='coordinate'   \
            level=0    \
            fixmate=1 \
            addupmarksupport=1 \
            tmpfile="${base_name}.tmp" \
            < "${base_name}.selected.bam" \
            > "${base_name}.sorted.bam"

        rm "${base_name}.selected.bam"

        # 15. Final conversion to CRAM
        scramble \
            -I bam \
            -O cram \
            -r "${reference_fasta}" \
            < "${base_name}.sorted.bam" \
            > "${output_cram_file}"

        rm "${base_name}.sorted.bam"
        """
}
