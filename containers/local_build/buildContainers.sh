# Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

echo "| --- buildContainers.sh -----------------------------------------------|"
echo "| USAGE: /bin/bash/ buildContainers.sh"
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

# echo "@ building amprecon_base_container.sif"
# singularity build --fakeroot amprecon_base_container.sif SingularityFile_amprecon_base_container

# echo "@ building amprecon_grc_tools.sif"
# singularity build --fakeroot amprecon_grc_tools.sif SingularityFile_amprecon_grc_tools

# echo "@ building amprecon_coi.sif"
# singularity build --fakeroot amprecon_coi.sif SingularityFile_amprecon_coi

echo "@ building amprecon_irods.sif"
singularity build --fakeroot amprecon_irods.sif SingularityFile_amprecon_irods

echo ":: DONE ::"
