clc
clearvars
close all

%heaps law.
load gene_strain_mat.mat
load strain_list.mat

N=250;
params=zeros(N,2);
yparams=params;
%X is for total, Y is for core
X=zeros(size(gene_strain_mat,2),N);
Y=X;
for I=1:N
    matrix=gene_strain_mat(:,randperm(size(gene_strain_mat,2)));
    
    y=true(size(matrix,1),1);
    x=false(size(matrix,1),1);
    C=zeros(size(matrix,2),1);
    D=C;
    for J=1:size(matrix,2)
        y=y&matrix(:,J);
        x=x|matrix(:,J);
        C(J)=sum(x);
        D(J)=sum(y);
    end
    P=polyfit(log([1:length(C)]),log(C),1);
    params(I,:)=P;
    X(:,I)=C;

    P=polyfit(log([1:length(D)]),log(D),1);
    yparams(I,:)=P;
    Y(:,I)=D;
    if mod(I,25)==0
        I
    end
end

C=flip(plasma(3));
C=C(1:2,:);
close all
f=figure;
x=1:size(gene_strain_mat,2);
plot(x,median(X,2),'Color',C(1,:),'LineWidth',2)
hold on
plot(x,prctile(X',5),'k--')
plot(x,prctile(X',95),'k--')

plot(x,median(Y,2),'Color',C(2,:),'LineWidth',2)
plot(x,prctile(Y',5),'k--')
plot(x,prctile(Y',95),'k--')


xlabel('Number of Strains')
ylabel('Number of Gene Clusters')
title('Heap''s Law Curve for the {\it B. fragilis} Pan Genome')
legend(['Median Cumlative Count (b = ',num2str(round(mean(params(:,1)),3)),')'], '','',...
    ['Median Conserved Count (b = ',num2str(round(mean(yparams(:,1)),3)),')'],'Location','east')

%compare to others
f=figure();

% doi.org/10.3389/fmicb.2019.00834 
%https://www.sciencedirect.com/science/article/pii/S0740002023001211
P=[exp(mean(params(:,2))),mean(params(:,1));904 .496; 1236 .329; 2404 .373; 5559 .375; 958 .435; 4022 .3491; 2772 .3395;  4300 .276; 2720 .301 ];
names={'Bacteroides fragilis','Streptococcus pneumoniea','Staphylococcus aureus','Salmonella entrica','Escherichia coli','Mycobacterium tuberculoisis','Pseudomonas aeruginosa','Acinetobacter baumanii','Bacillus subtilis','Lactiplantibacillus plantarum'};
[~,ind]=sort(P(:,1).*803.^P(:,2),'descend');
P=P(ind,:);
names=names(ind);
C=flip(plasma(length(P)));


hold on
for I=1:length(P)
    y2=P(I,1)*x.^P(I,2);
    if strcmp(names{I},'Bacteroides fragilis')
        plot(x,y2,'Color',C(I,:),'LineWidth',3)
        continue
    end
    plot(x,y2,'Color',C(I,:),'LineWidth',1)

end
legend(names, ...
    'Location','northwest')
title('Heap''s Law Curves')
xlabel('Number of Strains')
ylabel('Cumulative Pan-Genome Size')



%%
clc
clearvars
close all
%tsne to show no separation of MAGs and isolates
load gene_strain_mat.mat
load strain_list.mat



x=tsne(1*gene_strain_mat','Perplexity',5);

C=flip(plasma(11));


scatter(x(1:344,1),x(1:344,2),30,C(9,:),'filled')
hold on
scatter(x(345:end,1),x(345:end,2),30,C(3,:),'filled')

title('tSNE of Gene Presence')
xlabel('Dimension One (Unitless)')
ylabel('Dimension Two (Unitless)')
legend({'Isolates','MAGs'},'Location','northwest')


%%
%clustergram
clc
clearvars
close all
load gene_strain_mat.mat
% gene_strain_mat=gene_strain_mat(:,1:344);
d=(1-squareform(pdist(1*gene_strain_mat','jaccard')));
c=clustergram(d);
idx=cellfun(@str2num,c.ColumnLabels);
idy=cellfun(@str2num,c.RowLabels);
heatmap(d(idy,idx),'GridVisible','off','Colormap',flipud(plasma))
clim([0,1])

Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
title('Jaccard Similarity Between Strains')

xlabel('Strain')
ylabel('Strain')
