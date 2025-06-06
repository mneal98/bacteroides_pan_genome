clc
clearvars
close all

%correlate relative abundances to the ability to degrade certain n-linked
%glycans
data=importdata('strain_by_species_abundance.mat');
species=importdata('species_names.mat');
load rxn_strain_matrix.mat
load pan_model.mat
rxn_ind=find(strcmp(model.rxns,'G01419_decomposition'));
%repeat for
% G01419_decomposition 
% G02741_decomposition
% G03801_decomposition
% G11022_decomposition


[~,idx]=find(mean(data)>10^-5);
data=data(:,idx);
data(data==0)=min(data(data>0))/2;
species_names=species(idx);

%take the additive log-ratio relative to Bf
idx=find(contains(species_names,'Bacteroides fragilis'));
bfs_ab=sum(data(:,idx),2);
data(:,idx)=[];
species_names(idx)=[];
data=log2(data./bfs_ab);

strain_list=importdata('../strain_list.mat');
meta_G_strain_list=importdata("meta_G_strains.mat");
[~,idx]=intersect(strain_list,meta_G_strain_list);

gene_present=rxn_strain_matrix(rxn_ind,idx);
p_vals=zeros(size(species_names));
delta=p_vals;
for I=1:length(species_names)
    x=data(gene_present,I);
    y=data(~gene_present,I);
    [~,p]=ttest2(x,y);
    p_vals(I,1)=p;
    delta(I,1)=mean(x)-mean(y);
end
p_vals=mafdr(p_vals);
significant=find(p_vals<0.05);
p_vals=p_vals(significant);
delta=delta(significant);
species_names=species_names(significant);

for I=1:length(species_names)
    species_names{I}=extractBefore(species_names{I},' ');
end
[u, v]=groupcounts(species_names);
num_groups=sum(u>10);
C=flip(plasma(num_groups+3));
most_present_genera=v(u>10);

seen=false(size(species_names));
for I=1:num_groups
    ind=find(strcmp(species_names,most_present_genera{I}));
    seen(ind)=true;
end
plots=[];
plots(1)=scatter(delta(~seen),p_vals(~seen),75,C(2,:),'filled');
hold on
for I=1:num_groups
    ind=find(strcmp(species_names,most_present_genera{I}));
    plots(end+1)=scatter(delta(ind),p_vals(ind),75,C(I+2,:),'filled');
end
legend(flip(plots),flip(['Other';most_present_genera]),'Location','southwest')

title('Association Between Gut Species and Degradation Reaction Presence')
ylabel('P-Value')
xlabel('Log_2 Fold Change (With BfUbb / Without)')