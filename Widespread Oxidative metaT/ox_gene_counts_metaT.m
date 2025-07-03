clc
clearvars
close all

%here we sum up the CPM of the 10 identified genes that act as hallmarks of
%oxidative metabolism.
sample_gene_counts_CPM=importdata('sample_gene_counts_CPM.mat');
genes=importdata('gene_functions.mat');


labels={'Cytochrome Oxidase','Nitrite Reductase','NADH Dehydrogenase','ATP Synthase',...
        '6PGD','Citrate Synthase','Isocitrate Dehydrogenase',...
        'Succinyl-CoA Synthetase','Fumarase','Malate Dehydrogenase'};
%corresponding identifiers used in the data
gene_names={'Cytochrome bd','Nitroreductase','NADH-quinone','ATP synth',...
            '6-phosphogluconate dehydrogenase','Citrate synth','Isocitrate dehydrogenase',...
            'Succinate--CoA ligase','umarate hydratase','Malate dehydrog'
            };
ox_gene_counts=zeros(71,10);
for I=1:length(gene_names)
    ind=find(contains(genes,gene_names{I}));
    ox_gene_counts(:,I)=sum(sample_gene_counts_CPM(:,ind),2);
end

%percent of samples each gene is active in
percent_active=round(100*mean(ox_gene_counts>10))


C=flip(plasma(11));
[~,idx]=sort(mean(ox_gene_counts),'ascend');
ox_gene_counts=ox_gene_counts(:,idx);
labels=labels(idx);

b=boxplot(ox_gene_counts,'Labels',labels,'BoxStyle','filled','Colors',C,'Widths',0.5);


title('Expression of Respiratory Genes in Gut Metatranscriptomics')
ylabel('Expression level (CPM)')



%%
clc
clearvars
close all

%now we repeat, but compare to the CPM in the glucose-only sample
sample_gene_counts_CPM=importdata('sample_gene_counts_CPM.mat')+1;%add pseudocount
genes=importdata('gene_functions.mat');

glucose_only_data=table2cell(readtable('single_strain_translatomics_CPM.xlsx'));
glucose_only_cpm=cell2mat(glucose_only_data(:,2))+1;%add pseudocount
glucose_only_genes=glucose_only_data(:,1);

labels={'Cytochrome Oxidase','Nitrite Reductase','NADH Dehydrogenase','ATP Synthase',...
        '6PGD','Citrate Synthase','Isocitrate Dehydrogenase',...
        'Succinyl-CoA Synthetase','Fumarase','Malate Dehydrogenase'};
%corresponding identifiers used in the data
gene_names={'Cytochrome bd','Nitroreductase','NADH-quinone','ATP synth',...
            '6-phosphogluconate dehydrogenase','Citrate synth','Isocitrate dehydrogenase',...
            'Succinate--CoA ligase','umarate hydratase','Malate dehydrog'
            };
ox_gene_counts=zeros(71,10);
ox_gene_counts_isolate=zeros(10,1);
for I=1:length(gene_names)
    ind=find(contains(genes,gene_names{I}));
    ox_gene_counts(:,I)=sum(sample_gene_counts_CPM(:,ind),2);

    ox_gene_counts_isolate(I)=sum(glucose_only_cpm(contains(glucose_only_genes,gene_names{I})));
end

%percent of samples each gene is active in
percent_active=round(100*mean(ox_gene_counts>10))




C=flip(plasma(11));
expr_ratio=ox_gene_counts./(ox_gene_counts_isolate');
[~,idx]=sort(mean(expr_ratio),'ascend');

expr_ratio=expr_ratio(:,idx);
labels=labels(idx);
expr_ratio=log2(expr_ratio+1);

b=boxplot(expr_ratio,'Labels',labels,'BoxStyle','filled','Colors',C,'Widths',0.5);


title('Relative Expression of Respiratory Genes')
ylabel('Log_2 Fold Change (Gut / Glucose Alone)')



%%
clc
clearvars
close all

%make this data into PCA biplot
sample_gene_counts_CPM=importdata('sample_gene_counts_CPM.mat');
genes=importdata('gene_functions.mat');
labels={'Cytochrome Oxidase','Nitrite Reductase','NADH Dehydrogenase','ATP Synthase',...
        '6PGD','Citrate Synthase','Isocitrate Dehydrogenase',...
        'Succinyl-CoA Synthetase','Fumarase','Malate Dehydrogenase'};

gene_names={'Cytochrome bd','Nitroreductase','NADH-quinone','ATP synth',...
            '6-phosphogluconate dehydrogenase','Citrate synth','Isocitrate dehydrogenase',...
            'Succinate--CoA ligase','umarate hydratase','Malate dehydrog'
            };
ox_gene_counts=zeros(71,10);
for I=1:length(gene_names)
    ind=(contains(genes,gene_names{I}));
    ox_gene_counts(:,I)=sum(sample_gene_counts_CPM(:,ind),2);
end

a=sum(ox_gene_counts,2);
[~,x,~,~,e]=pca(sample_gene_counts_CPM);


x=flip(x,1);
a=flip(a);
scatter(x(:,1),x(:,2),50,a,'filled');
colormap(flipud(plasma))
[~,ind]=max(a);
hold on
scatter(x(ind,1),x(ind,2),50,a(ind),'filled');
[~,ind]=min(a);
scatter(x(ind,1),x(ind,2),50,a(ind),'filled');


xlabel(['Component 1 (',num2str(round(e(1))),'%)'])
ylabel(['Component 2 (',num2str(round(e(2))),'%)'])

title('PCA of {\itin vivo Bacteroides fragilis} Expression')
legend({'','High Respiratory Gene Activity','Low Respiratory Gene Activity'})

%%
%and a heatmap
clc
clearvars
close all

sample_gene_counts_CPM=importdata('sample_gene_counts_CPM.mat')+1;
genes=importdata('gene_functions.mat');
labels={'Cytochrome Oxidase','Nitrite Reductase','NADH Dehydrogenase','ATP Synthase',...
        '6PGD','Citrate Synthase','Isocitrate Dehydrogenase',...
        'Succinyl-CoA Synthetase','Fumarase','Malate Dehydrogenase'};

gene_names={'Cytochrome bd','Nitroreductase','NADH-quinone','ATP synth',...
            '6-phosphogluconate dehydrogenase','Citrate synth','Isocitrate dehydrogenase',...
            'Succinate--CoA ligase','umarate hydratase','Malate dehydrog'
            };
ox_gene_counts=zeros(71,10);
for I=1:length(gene_names)
    ind=(contains(genes,gene_names{I}));
    ox_gene_counts(:,I)=sum(sample_gene_counts_CPM(:,ind),2);
end
ox_gene_counts=log10(ox_gene_counts);

[~,idx]=sort(sum(ox_gene_counts,2));
[~,idy]=sort(sum(ox_gene_counts,1));
% idx=1:71;
% idy=1:10;
heatmap(ox_gene_counts(idx,idy),'GridVisible','off','Colormap',flipud(plasma),...
    'XData',labels(idy))
Ax = gca;
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
ylabel('Gut Metatranscriptomic Sample')
title('Expression of Respiratory Genes log_1_0(CPM)')