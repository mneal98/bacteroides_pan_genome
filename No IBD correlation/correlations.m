clc
clearvars
close all

%In this script, we try to separate IBD from healthy associated strains.
load gene_strain_mat.mat
load rxn_strain_matrix.mat
load strain_list.mat
load pan_model.mat
C=importdata('../Diversity/all_growth_C.mat');
N=importdata('../Diversity/all_growth_N.mat');
all_growth=[C,N]';

opts = detectImportOptions('../Strain Metadata.xlsx');
opts=setvartype(opts,'char');
metadata=table2cell(readtable('../Strain Metadata.xlsx',opts));
healthy=metadata(strcmp(metadata(:,2),'Healthy'),1);
[~, healthy, ~]=intersect(strain_list,healthy);
IBD=metadata(strcmp(metadata(:,2),'IBD'),1);
[~, IBD, ~]=intersect(strain_list,IBD);
infection=metadata(strcmp(metadata(:,2),'Infection'),1);
[~, infection, ~]=intersect(strain_list,infection);

labels=[repmat({'Healthy'},size(healthy));repmat({'IBD'},size(IBD));repmat({'Infection'},size(infection))];
%binary growth data
data=1*([all_growth(:,healthy),all_growth(:,IBD),all_growth(:,infection)]>0.01);

num_mets = size(data, 1);
pvals = zeros(num_mets, 1);
groups = unique(labels); 
for i = 1:num_mets
    met = data(i, :);  % metabolite utilization across strains
    % Build 2x3 contingency table: rows = [present; absent], cols = [group 1; 2; 3]
    table = zeros(2, length(groups));
    for g = 1:length(groups)
        idx = find(strcmp(labels,groups{g}));
        table(1, g) = sum(met(idx) == 1);  % present
        table(2, g) = sum(met(idx) == 0);  % absent
    end
    
    % Chi-squared test
    [~,~,p]=chi2ind(table,0.0);
    
    pvals(i) = p;
end
pvals(isnan(pvals))=1;
% Multiple testing correction (FDR)
[~, adj_pvals] = mafdr(pvals);
% Get significant indecies
sig_idx=find(adj_pvals < 0.05);
'Number of distinguishing C/N sources'
length(sig_idx)








%and repeat for reactions
data=1*[rxn_strain_matrix(:,healthy),rxn_strain_matrix(:,IBD),rxn_strain_matrix(:,infection)];
num_rxns = size(data, 1);
pvals = zeros(num_rxns, 1);
groups = unique(labels); 
for i = 1:num_rxns
    met = data(i, :);  % metabolite utilization across strains
    % Build 2x3 contingency table: rows = [present; absent], cols = [group 1; 2; 3]
    table = zeros(2, length(groups));
    for g = 1:length(groups)
        idx = find(strcmp(labels,groups{g}));
        table(1, g) = sum(met(idx) == 1);  % present
        table(2, g) = sum(met(idx) == 0);  % absent
    end
    % Chi-squared test
    [~,~,p]=chi2ind(table,0.0);
    pvals(i) = p;
end
pvals(isnan(pvals))=1;
% Multiple testing correction (FDR)
[~, adj_pvals] = mafdr(pvals);
% Get significant indecies
sig_idx=find(adj_pvals < 0.05);
'Significant Reactions'
[model.rxns(sig_idx),model.rxnNames(sig_idx)]