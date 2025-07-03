clc
clearvars
close all

%correlate relative abundances to the presence/abscence of certain genes.
%first lets look at the secreted ubiquitin like protein
data=importdata('strain_by_species_abundance.mat');
species=importdata('species_names.mat');
pan_genome_sequences=importdata('pan_genome_sequences.mat');
load gene_strain_mat.mat

gene_ind=35582; %where ubiquitin is found
pan_genome_sequences(gene_ind).Header

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

gene_present=gene_strain_mat(gene_ind,idx);
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
species_names_orig=species_names;
for I=1:length(species_names)
    species_names{I}=extractBefore(species_names{I},' ');
end
[u, v]=groupcounts(species_names);
num_groups=sum(u>1);
C=flip(plasma(num_groups+3));
most_present_genera=v(u>1);

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
for I=1:length(most_present_genera)
    most_present_genera{I}=['{\it',most_present_genera{I},'}'];
end
legend(flip(plots),flip(['Other';most_present_genera]),'Location','southwest')

title('Association Between Gut Species and BfUbb Presence')
ylabel('P-Value')
xlabel('Log_2 Fold Change (With BfUbb / Without)')



%%
clc
clearvars
close all

%repeat for the bacteroides fragilis toxin, bft
data=importdata('strain_by_species_abundance.mat');
species=importdata('species_names.mat');
pan_genome_sequences=importdata('pan_genome_sequences.mat');
load gene_strain_mat.mat

gene_ind=9903;

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

gene_present=gene_strain_mat(gene_ind,idx);
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
species_names_orig=species_names;
for I=1:length(species_names)
    species_names{I}=extractBefore(species_names{I},' ');
    if contains(species_names{I},'[')
        a=extractBetween(species_names{I},'[',']');
        species_names{I}=a{1};
    end
end


num_groups=3;
C=flip(plasma(num_groups+5));
[u, v]=groupcounts(species_names);
[a, b]=maxk(u,num_groups);
most_present_genera=sort(v(b));

[u, v]=groupcounts(species_names(delta>0));
[a, b]=maxk(u,2);
most_present_positive_genera=sort(v(b));

seen=false(size(species_names));
for I=1:num_groups
    ind=find(strcmp(species_names,most_present_genera{I}));
    seen(ind)=true;
end
for I=1:2
    ind=find(strcmp(species_names,most_present_positive_genera{I}));
    seen(ind)=true;
end

plots=[];
plots(1)=scatter(delta(~seen),p_vals(~seen),75,C(2,:),'filled');
hold on
for I=1:num_groups
    ind=find(strcmp(species_names,most_present_genera{I}));
    plots(end+1)=scatter(delta(ind),p_vals(ind),75,C(I+2,:),'filled');
end
for I=1:2
    ind=find(strcmp(species_names,most_present_positive_genera{I}));
    plots(end+1)=scatter(delta(ind),p_vals(ind),75,C(I+2+num_groups,:),'filled');
end
for I=1:length(most_present_genera)
    most_present_genera{I}=['{\it',most_present_genera{I},'}'];
end
for I=1:length(most_present_positive_genera)
    most_present_positive_genera{I}=['{\it',most_present_positive_genera{I},'}'];
end
legend(flip(plots),flip(['Other';most_present_genera;most_present_positive_genera]),'Location','northwest')
title('Association Between Gut Species and BFT Toxin Presence')
ylabel('P-Value')
xlabel('Log_2 Fold Change (With BFT / Without)')






%%
clc
clearvars
close all

%correlate relative abundances to the presence/abscence of certain genes.
%first lets look at the secreted ubiquitin like protein
data=importdata('strain_by_species_abundance.mat');
species=importdata('species_names.mat');
pan_genome_sequences=importdata('pan_genome_sequences.mat');
load gene_strain_mat.mat

gene_ind=35582; %where ubiquitin is found
pan_genome_sequences(gene_ind).Header

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

gene_present=gene_strain_mat(gene_ind,idx);
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
num_groups=sum(u>1);
C=flip(plasma(num_groups+3));
most_present_genera=v(u>1);

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

title('Association Between Gut Species and BfUbb Presence')
ylabel('P-Value')
xlabel('Log_2 Fold Change (With BfUbb / Without)')



