The folders contain scripts that make the figures for each main section.

For utilizing the models, the main files are:
1) pan_model.mat. This is the pan-genome model in matlab format. The reactions and metabolites are also listed in the spreadsheet pan_model.xlsx.
2) strain_list.mat / Strain Identifiers.xlsx. This lists the strains used by their BVBRC or MGNIFY Identifiers.
3) rxn_strain_matrix.mat. This binary matrix encodes which strains have which reactions. Element i,j is 1 if reaction i from the the pan_model is present in strain j from strain_list.
4) pan_genome_sequences.mat. This lists all the genes in the pan-genome, with identifiers and amino acid sequences.
5) gene_strain_mat.mat. This binary matrix encodes which strains have which genes from pan_genome_sequences.

With these files, you may reconstruct a strain specific model by taking the pan_model and deleting the reactions not present in that strain, e.g. by using the command:
strain_model = removeRxns( pan_model, pan_model.rxns(rxn_strain_matrix(:,I)==0) ).

The genes of the pan_model are written in the form "x(index)" where "index" refers to an entry of the pan_genome_sequences. For example, Thymidylate synthase (TMDS) has the GPR x(18998), which indicates that pan_genome_sequences gene number 18998 is needed to perform this reaction. The rxn_strain_matrix encodes the pre-computed boolean evaluations of these GPRs using the genetic inventory of each strain. Additionally, "GAP" indicates necessary reactions with no genes found, "SPONT" indicated spontaneous reactions, and "EX" indicated exchange reactions.


The processed metagenomic data is found in the file 'strain_by_species_abundance.mat', with the associated B. fragilis strains in the file 'meta_G_strains.mat' and the co-occuring species listed in species_names.mat.

The processed transcriptomic data from the in vitro experiments is found in 'single_strain_translatomics_CPM.xlsx,' and the metaT is found in 'sample_gene_counts_CPM.mat,' 'gene_functions.mat,' and metadata in 'disease_status.mat' and 'metatranscriptomic_metadata.xlsx.'
