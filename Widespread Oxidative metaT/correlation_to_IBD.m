clc
clearvars
close all

%here we identify any associations between IBD/nonIBD staus and metaT
%results

disease_status=importdata('disease_status.mat');
sample_gene_counts_CPM=importdata('sample_gene_counts_CPM.mat');
gene_functions=importdata('gene_functions.mat');

IBD=sample_gene_counts_CPM(~strcmp(disease_status,'nonIBD'),:);
nonIBD=sample_gene_counts_CPM(strcmp(disease_status,'nonIBD'),:);

[~,pvals]=ttest2(IBD,nonIBD);
%FDR adjust
FDR=mafdr(pvals);
ind=find(FDR<.01);

%here are the few significant genes. 
significant_genes=gene_functions(ind)