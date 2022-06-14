
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

echo "@ building biobambam2.sif"
singularity build --fakeroot biobambam2.sif SingularityFile_biobambam2

echo "@ building staden.sif"
singularity build --fakeroot staden.sif SingularityFile_staden

echo "@ pulling samtools_1.15.sif"
singularity pull docker://staphb/samtools:1.15

echo "@ pulling bwa_0.7.17.sif"
singularity pull docker://staphb/bwa:0.7.17

echo ":: DONE ::"
