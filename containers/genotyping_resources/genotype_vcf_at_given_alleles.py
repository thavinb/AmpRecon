#!/usr/bin/env python3

import sys
from optparse import OptionParser
import errno

import os
import numpy as np
import copy
import vcf
import allel
import datetime
import collections
from subprocess import call

__version__ = "0.0.3"

if __name__ == '__main__':

    usage = 'usage: %prog [options]'
    description = "Genotype a vcf at a specific set of alleles"
    epilog = """
This genotypes a vcf (typically the output of running GenotypeGVCFs with -allSites) at a specific set of alleles from a second
vcf file.

Examples:

    genotype_vcf_at_given_alleles_20180608.py \
    --gvcf_fn /lustre/scratch116/malaria/pfalciparum/output/3/b/5/7/338034/1_gatk_genotype_gvcfs_gatk3_v2/gatk.genotyped.vcf.gz \
    --alleles_fn /lustre/scratch118/malaria/team112/pipelines/resources/pfalciparum/pf_6x_genotyping.vcf.gz \
    --genotyped_fn QV0090-C.pf_6x.vcf \
    --GenotypeGVCFs_format {{sample}}.GenotypeGVCFs.{{alleles}}.vcf.gz \
    --java_command /software/jre1.8.0_131/bin/java -Xms5625m -server -XX:+UseSerialGC \
    --gatk_jar /nfs/team112_internal/production/tools/bin/gatk/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar \
    --bcftools /nfs/users/nfs_r/rp7/miniconda3/envs/biipy/bin/bcftools \
    --bgzip /nfs/team112_internal/production/tools/bin/bgzip \
    --verbose

Version: {version}
""".format(version=__version__)

    OptionParser.format_epilog = lambda self, formatter: self.epilog
    parser = OptionParser(usage=usage, description=description, epilog=epilog)
    parser.add_option('-i', '--gvcf_fn', dest='gvcf_fn', help='input gvcf filename [/lustre/scratch116/malaria/pfalciparum/output/c/0/3/1/65458/1_gatk_haplotype_caller_gatk3_v3/gatk.vcf.gz]', default='/lustre/scratch116/malaria/pfalciparum/output/c/0/3/1/65458/1_gatk_haplotype_caller_gatk3_v3/gatk.vcf.gz')
    parser.add_option('-a', '--alleles_fn', dest='alleles_fn', help='alleles filename [/lustre/scratch118/malaria/team112/pipelines/resources/pfalciparum/pf_6x_genotyping.vcf.gz]', default='/lustre/scratch118/malaria/team112/pipelines/resources/pfalciparum/pf_6x_genotyping.vcf.gz')
    parser.add_option('-o', '--genotyped_fn', dest='genotyped_fn', help='output genotyped filename [QV0090-C.pf_6x.vcf]', default='QV0090-C.pf_6x.vcf')
    parser.add_option('-f', '--GenotypeGVCFs_format', dest='GenotypeGVCFs_format', help='format for interim output of GenotypeGVCFs [{{sample}}.GenotypeGVCFs.{{alleles}}.vcf.gz]', default='{{sample}}.GenotypeGVCFs.{{alleles}}.vcf.gz')
    parser.add_option('-j', '--java_command', dest='java_command', help='java command', default='/software/jre1.8.0_131/bin/java -server -XX:+UseSerialGC')
    parser.add_option('-m', '--java_memory', dest='java_memory', help='java momory', default='5625m')
    parser.add_option('-r', '--ref_fasta', dest='ref_fasta', help='Reference fasta file [/lustre/scratch116/malaria/pfalciparum/resources/Pfalciparum.genome.fasta]', default='/lustre/scratch116/malaria/pfalciparum/resources/Pfalciparum.genome.fasta')
    parser.add_option('-g', '--gatk_jar', dest='gatk_jar', help='GATK jar file', default='/nfs/team112_internal/production/tools/bin/gatk/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar')
    parser.add_option('-b', '--bcftools', dest='bcftools', help='bcftools executable', default='/nfs/users/nfs_r/rp7/miniconda3/envs/biipy/bin/bcftools')
    parser.add_option('-z', '--bgzip', dest='bgzip', help='bgzip executable', default='/nfs/team112_internal/production/tools/bin/bgzip')
    parser.add_option('-x', '--remove_genotype_gvcf_file', action='store_true', dest='remove_genotype_gvcf_file', default=False)
    parser.add_option('-w', '--overwrite', action='store_true', dest='overwrite', default=False)
    parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=True)
    options, args = parser.parse_args()
    options.java_command = f"{options.java_command} -Xmx{options.java_memory} -Xms{options.java_memory}"

    try:

        def genotype_vcf_at_given_alleles(
            gvcf_fn = '/lustre/scratch116/malaria/pfalciparum/output/c/0/3/1/65458/1_gatk_haplotype_caller_gatk3_v3/gatk.vcf.gz', # QV0090-C
            alleles_fn='/lustre/scratch118/malaria/team112/pipelines/resources/pfalciparum/pf_6x_genotyping.vcf.gz',
            genotyped_fn='QV0090-C.pf_6x.vcf',
            GenotypeGVCFs_format='{sample}.GenotypeGVCFs.{alleles}.vcf.gz',
            ref_fasta='/lustre/scratch116/malaria/pfalciparum/resources/Pfalciparum.genome.fasta',
            java_command='/software/jre1.8.0_131/bin/java -Xmx5625m -Xms5625m -server -XX:+UseSerialGC',
            gatk_jar='/nfs/team112_internal/production/tools/bin/gatk/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar',
            bcftools='/nfs/users/nfs_r/rp7/miniconda3/envs/biipy/bin/bcftools',
            bgzip='/nfs/team112_internal/production/tools/bin/bgzip',
            remove_genotype_gvcf_file=False,
            overwrite=False,
            verbose=True
        ):
            if not os.path.exists("%s.gz" % genotyped_fn) or overwrite:
                sample = vcf.Reader(filename=gvcf_fn).samples[0]
                if(verbose):
                    print(sample)

                GenotypeGVCFs_fn = GenotypeGVCFs_format.format(sample=sample, alleles=os.path.basename(alleles_fn).split('.')[0])
                if(verbose):
                    print(gvcf_fn)
                    print(GenotypeGVCFs_fn)
                    print()
                if not os.path.exists(GenotypeGVCFs_fn) or overwrite:
                    call(
                        java_command.split(' ') + [
                            "-jar",  gatk_jar, "-T", "GenotypeGVCFs", "-R", ref_fasta, "--variant", gvcf_fn,
                            "-allSites", "-L", alleles_fn, "-o", GenotypeGVCFs_fn
                        ]
                    )
                else:
                    print("File %s already exists - using this. If you would like recreate, rerun with --overwrite" % (GenotypeGVCFs_fn))

                all_sites_reader = vcf.Reader(filename=GenotypeGVCFs_fn)
                alleles_reader = vcf.Reader(filename=alleles_fn)
                vcf_writer = vcf.Writer(open(genotyped_fn, 'w'), all_sites_reader)
                current_chrom = ''
                same_position = True
                for alleles_record in alleles_reader:
                    if verbose:
                        if alleles_record.CHROM != current_chrom:
                            current_chrom = alleles_record.CHROM
                            print("Processing %s" % current_chrom)
                    if same_position:
                        all_sites_record = next(all_sites_reader)
                    if (all_sites_record.CHROM == alleles_record.CHROM) and (all_sites_record.POS == alleles_record.POS):
                        same_position = True
                        # Same alleles
                        if (all_sites_record.REF == alleles_record.REF) and (all_sites_record.ALT == alleles_record.ALT):
                            vcf_writer.write_record(all_sites_record)
                        # Different alleles
                        # Here we go through each alt allele in the alleles vcf, and ensure we have the correct data for that allele in output
                        else:
                            # See comment from Dawe at https://github.com/jamescasbon/PyVCF/issues/82
                            f_keys = all_sites_record.samples[0].data._fields
                            f_vals = [all_sites_record.samples[0].data[vx] for vx in range(len(all_sites_record.samples[0].data))]
                            handy_dict = dict(zip(f_keys, f_vals))

                            output_record = copy.deepcopy(all_sites_record)
                            output_record.ALT = alleles_record.ALT
                            output_record.REF = alleles_record.REF
                            number_of_alts = len(output_record.ALT)
                            # Some rows have DP but not AD. Where AD is missing we will use DP
                            if 'AD' in handy_dict.keys():
                                original_AD = copy.deepcopy(handy_dict['AD'])
                            else:
                                if 'DP' in handy_dict.keys():
                                    original_AD = [copy.deepcopy(handy_dict['DP'])]
                                else:
                                    original_AD = [0]
                                temp = list(f_keys)
                                temp.insert(1, 'AD')
                                f_keys = tuple(temp)
                                orig_format = output_record.FORMAT.split(':')
                                orig_format.insert(1, 'AD')
                                output_record.FORMAT = ':'.join(orig_format)
                            handy_dict['AD'] = [original_AD[0]] + ([0] * (number_of_alts))

                            output_record.samples[0].data = collections.namedtuple('CallData', f_keys)

                            for i, alleles_alt in enumerate(alleles_record.ALT):
                                if alleles_alt in all_sites_record.ALT:
                                    if verbose:
                                        print("Different alleles at: ", all_sites_record, alleles_record)
                                    j = np.where(np.array(all_sites_record.ALT) == alleles_alt)[0][0]
                                    GT_alleles_encodoing = str(i+1)
                                    GT_all_sites_encodoing = str(j+1)
                                    handy_dict['GT'] = handy_dict['GT'].replace(GT_all_sites_encodoing, GT_alleles_encodoing)
                                    handy_dict['AD'][i+1] = original_AD[j+1]
                            # If any genotypes for allele not in alleles vcf, set genotype to missing
                            parts=handy_dict['GT'].split('/')
                            if len(parts)==1: # Sometimes the separator is '|'... (could use re.split() but why bother)
                                parts = handy_dict['GT'].split('|')
                            temp_GTs = np.array(parts)
                            if not np.any(temp_GTs == '.') and np.any(temp_GTs.astype(int) > number_of_alts):
                                if verbose:
                                    print("Genotype for non-existent allele at: ", all_sites_record, alleles_record, handy_dict['GT'], ". Setting to missing.")
                                handy_dict['GT'] = './.'
                            # If too many PLs after removing alleles, reduce to correct PLs
                            expected_number_of_PLs = int((number_of_alts+1)*(number_of_alts+2)/2)
                            if 'PL' in handy_dict and len(handy_dict['PL']) > expected_number_of_PLs:
                                if verbose:
                                    print("Too many PLs at: ", all_sites_record, alleles_record, handy_dict['GT'], ". Setting to", expected_number_of_PLs)
                                handy_dict['PL'] = handy_dict['PL'][0:expected_number_of_PLs]
                            if 'PL' in handy_dict and len(handy_dict['PL']) < expected_number_of_PLs:
                                if verbose:
                                    print("Not enough PLs at: ", all_sites_record, alleles_record, handy_dict['GT'], ". Increasing to", expected_number_of_PLs)
                                max_PL = np.max(handy_dict['PL'])
                                extra_PLs = [max_PL] * (expected_number_of_PLs - len(handy_dict['PL']))
                                handy_dict['PL'].extend(extra_PLs)
                                # handy_dict.pop('PL')
                                # orig_format = output_record.FORMAT.split(':')
                                # orig_format.remove('PL')
                                # output_record.FORMAT = ':'.join(orig_format)
                                # temp = list(f_keys)
                                # temp.remove('PL')
                                # f_keys = tuple(temp)
                            new_vals = [handy_dict[x] for x in f_keys]
                            # finally set CallData
                            output_record.samples[0].data = output_record.samples[0].data._make(new_vals)

                            vcf_writer.write_record(output_record)
                    elif (all_sites_record.CHROM > alleles_record.CHROM) or (all_sites_record.POS > alleles_record.POS):
                        same_position = False
                    # Record missing in sample file, and hence enter information from alleles file and give an empty genotype
                        output_record = copy.deepcopy(alleles_record)
                        number_of_alts = len(output_record.ALT)
                        # We want to remove all INFO fields in cases where there is no record in output of GenotypeGCVFs
                        output_record.INFO = {}
                        output_record.FILTER = '.'
                        output_record.FORMAT = 'GT:AD:DP'
                        CallData = collections.namedtuple('CallData', ['GT', 'AD', 'DP'])
                        output_record.samples = [vcf.model._Call(output_record, sample, CallData(GT='./.', AD=([0] * (number_of_alts + 1)), DP=0))]
                        vcf_writer.write_record(output_record)
                    elif (all_sites_record.CHROM < alleles_record.CHROM) or (all_sites_record.POS < alleles_record.POS):
                        raise ValueError('There is a record at position %s:%d in genotype file that is not in alleles file. This should not happen!' % (all_sites_record.CHROM, all_sites_record.POS))
                    else:
                        raise Exception('Record in genotypes file neither before or after record in alleles file. This makes no sense. Maybe your files are in a weird order?')

                vcf_writer.close()

                # bgzip and index vcf
                call([bgzip, "-f", genotyped_fn])
                call([bcftools, "index", "--tbi", "%s.gz" % genotyped_fn])

                if remove_genotype_gvcf_file:
                    call(["rm", GenotypeGVCFs_fn])
                    call(["rm", "%s.tbi" % GenotypeGVCFs_fn])

                print("done")
            else:
                print("File %s.gz already exists - abandoning. If you would like recreate, rerun with --overwrite" % (genotyped_fn))


        genotype_vcf_at_given_alleles(
            gvcf_fn=options.gvcf_fn,
            alleles_fn=options.alleles_fn,
            genotyped_fn=options.genotyped_fn,
            GenotypeGVCFs_format=options.GenotypeGVCFs_format,
            ref_fasta=options.ref_fasta,
            java_command=options.java_command,
            gatk_jar=options.gatk_jar,
            bcftools=options.bcftools,
            bgzip=options.bgzip,
            remove_genotype_gvcf_file=options.remove_genotype_gvcf_file,
            overwrite=options.overwrite,
            verbose=options.verbose,
        )

    except IOError as e:
        if e.errno == errno.EPIPE:
            pass # ignore broken pipe
        else:
            raise
