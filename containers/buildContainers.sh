
echo "| --- setupContainers.sh -----------------------------------------------|"
echo "| USAGE: /bin/bash/ setupContainers.sh"
echo "| This scripts setups containers to be used on ampseq-pipeline "
echo "| Be sure all required singularity recipies are present at working dir."
echo "| ----------------------------------------------------------------------|"

# --- SETUP ERROR HANDLING ----------------------------------------------------
# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

echo "@ building bambi.sif"
singularity build --fakeroot bambi.sif SingularityFile_bambi

echo "@ building biobambam2-2.0.79.sif"
singularity build --fakeroot biobambam2-2.0.79.sif biobambam2-2.0.79.def

echo "@ building staden.sif"
singularity build --fakeroot staden.sif SingularityFile_staden

echo "@ building mergeheaders.sif"
singularity build --fakeroot mergeHeaders.sif SingularityFile_mergeHeaders

echo "@ pulling samtools_1.15.sif"
singularity pull docker://staphb/samtools:1.15

echo "@ pulling samtools_1.9.sif" # same of current production pipeline
singularity pull docker://staphb/samtools:1.9


echo "@ building picard.sif"
singularity build --fakeroot picard.sif SingularityFile_picard

echo "@ building bwa.sif"
singularity build --fakeroot bwa.sif SingularityFile_bwa

echo "@ building pythonBox.sif"
singularity build --fakeroot pythonBox.sif SingularityFile_pythonBox

echo "@ building genotyping.sif"
singularity build --fakeroot genotyping.sif SingularityFile_genotyping

echo "@ building bcftools.sif"
singularity build --fakeroot bcftools.sif SingularityFile_bcftools

echo "@ building grc_tools.sif"
singularity build --fakeroot grc_tools.sif SingularityFile_grctools

echo "@ building coi.sif"
singularity build --fakeroot coi.sif SingularityFile_coi

echo "@ building pyvcf.sif"
singularity build --fakeroot pyvcf.sif SingularityFile_pyvcf

echo ":: DONE ::"
