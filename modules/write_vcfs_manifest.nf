process write_vcfs_manifest {

    input:
        // this script assumes "all at once" will be provided to this process 
        val(IDs_list)
        val(vcf_paths_list)

    output:
        file("${mnf_out_nm}")

    script:
    mnf_out_nm = "lanelet_vcf_manifest.csv"
    """
#!/usr/bin/python3

import csv

def write_csv(id_list, vcf_path_list, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["ID", "vcf_path"])
        for id, vcf_path in zip(id_list, vcf_path_list):
            writer.writerow([id, vcf_path])

    #write_csv(id_list, vcf_path_list, output_file)
    # get lists
    # id_list = ["id_1", "id_2"]
id_list = list("${IDs_list}".strip("[]").replace(" ", "").split(","))

#vcf_path_list = ["vcf_path_1", "vcf_path_2"]
vcf_path_list =  list("${vcf_paths_list}".strip("[]").replace(" ", "").split(","))

write_csv(id_list, vcf_path_list, "${mnf_out_nm}")
    """
}
