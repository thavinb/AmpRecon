
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

echo "@ building default.sif"
singularity build --fakeroot default.sif SingularityFile_default

echo "@ building genotyping.sif"
singularity build --fakeroot genotyping.sif SingularityFile_genotyping

echo "@ building grc_tools.sif"
singularity build --fakeroot grc_tools.sif SingularityFile_grctools

echo "@ building coi.sif"
singularity build --fakeroot coi.sif SingularityFile_coi

echo ":: DONE ::"
