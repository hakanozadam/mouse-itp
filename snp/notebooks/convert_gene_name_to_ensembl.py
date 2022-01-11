name_table_file = "./name_correspondences.tsv"

def get_gene_to_dict(name_table_file, name_type = "transcript"):
    if name_type == "transcript":
        pick_index = 2
    else:
        pick_index = 1

    result = dict()

    with open(name_table_file, "rt") as input_stream:
        for line in input_stream:
            contents = line.split()

            if len(contents) < 3:
                continue

            result[contents[0]] = contents[pick_index]

    return result


gene_to_ensembl_transcript_dict = get_gene_to_dict(name_table_file, name_type = "transcript")
gene_to_ensembl_gene_dict       = get_gene_to_dict(name_table_file, name_type = "gene")

# Returns / prints the Ensembl gene name
print(gene_to_ensembl_gene_dict["Nin"] )

# returns the Ensembl transcript name
print(gene_to_ensembl_transcript_dict["Nin"])
