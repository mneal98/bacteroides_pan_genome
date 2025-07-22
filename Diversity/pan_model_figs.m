clc
clearvars
close all
load rxn_strain_matrix.mat
load gene_strain_mat.mat
load pan_model.mat

%bar chart of the genes rxns mets

r2=rxn_strain_matrix(~startsWith(model.rxns,'EX_'),:);
R=sum(r2);
mean_rxns=round(mean(R))

x=false(length(gene_strain_mat),1);
for I=1:length(model.rxns)
    gpr=model.grRules{I};
    inds=extractBetween(gpr,'x(',')');
    inds=cellfun(@str2num,inds);
    x(inds)=1;
end
gsm2=gene_strain_mat.*x;
G=sum(gsm2);
mean_genes=round(mean(G))

for I=1:700
    m2=removeRxns(model,model.rxns(~rxn_strain_matrix(:,I)));
    M(I)=length(m2.mets);
end
mean_mets=round(mean(M))

X=[G;R;M]';

C=flip(plasma(4));
labels={'Genes','Reactions','Metabolites'};
b=violinplot(X,labels);
for I=1:3
    b(I).ViolinPlot.FaceColor=C(I,:);
    b(I).ViolinPlot.EdgeColor=C(I,:);
    b(I).ScatterPlot.MarkerFaceColor=C(I,:);
    b(I).ScatterPlot.MarkerFaceAlpha=0.0;
end

title('Pan-Reactome Content Distribution')
ylabel('Count')
set(gca,'FontName','Arial')
set(gca,'FontSize',14)


%%
clc
clearvars
close all
%core accessory by subsystem
load pan_model.mat
load rxn_strain_matrix.mat

core=find(sum(rxn_strain_matrix')>=(700*.99));
accessory=find(sum(rxn_strain_matrix')<(700*.99));


lipids={'Cell Envelope Biosynthesis','Membrane Lipid Metabolism','Fatty Acid Metabolism','Glycerophospholipid Metabolism','Mycolic acid pathway','Pantothenate and CoA biosynthesis','Steroid hormone biosynthesis'};
amino={'Alanine, aspartate and glutamate metabolism','Amino Acid Metabolism','Arginine and Proline Metabolism','Cysteine and methionine metabolism','Glycine, Serine, and Threonine Metabolism','Histidine Metabolism','Peptidoglycan Metabolism','Peptidoglycan biosynthesis','Phenylalanine, tyrosine and tryptophan biosynthesis','Selenoamino acid metabolism','Threonine and Lysine Metabolism','Valine, Leucine, and Isoleucine Metabolism','tRNA Charging','Glutathione metabolism','D-Glutamine and D-glutamate metabolism'};
carbs={'Amino sugar and nucleotide sugar metabolism','Alternate Carbon Metabolism','Citric Acid Cycle','Central Metabolism','Fructose and mannose metabolism','Galactose metabolism','Glycolysis / Gluconeogenesis','Pentose Phosphate Pathway','Starch and sucrose metabolism','Pyruvate Metabolism'};
transport={'Transport'};
secondary={'Biosynthesis of secondary metabolites','Cofactor and Prosthetic Group Biosynthesis','Drug metabolism - other enzymes','Folate Metabolism','Nicotinate and nicotinamide metabolism','Porphyrin and chlorophyll metabolism','Taurine and hypotaurine metabolism','Terpenoid backbone biosynthesis','Ubiquinone biosynthesis'};
nucleic={'Nucleotide Salvage Pathway','Purine and Pyrimidine Biosynthesis','Purine metabolism','Pyrimidine Metabolism','Pyrimidine metabolism'};
groups={lipids,amino,carbs,transport,secondary,nucleic};
labels={'Lipids','Amino Acids','Carbohydrates','Transport','Secondary Metabolites','Nucleic Acids'};

%core, accessory
counts=zeros(length(groups),2);
for I=1:length(groups)
    lst=groups{I};
    inds=[];
    for J=1:length(lst)
        inds=[inds;find(strcmp(model.subSystems,lst{J}))];
    end
    counts(I,1)=length(intersect(core,inds));
    counts(I,2)=length(intersect(accessory,inds));
end
[~,ind]=sort(sum(counts,2),'descend');
counts=counts(ind,:);
labels=labels(ind);

C=flip(plasma(4));
X=categorical(labels);
X=reordercats(X,labels);

b=bar(X,counts,'stacked','FaceColor','flat')
for I=1:length(b)
    b(I).CData=0* b(I).CData+C(I,:);
end
legend({'Core ','Accessory'},'Location','northeast')
ylabel('Count')
title('Core and Accessory Reactions by Category')
set(gca,'FontName','Arial')
set(gca,'FontSize',14)

%%
clc
clearvars
close all
%core accessory growth producing metabolites
C=importdata("all_growth_C.mat");
N=importdata("all_growth_N.mat");
load pan_model.mat

exchanges=model.rxns(startsWith(model.rxns,'EX_'));

can_grow=sum(C>0.001);
core_C=sum(can_grow>=(700*.99));
accessory_C=sum((can_grow<(700*.99)).*(can_grow>0));

can_grow=sum(N>0.001);
core_N=sum(can_grow>=(700*.99));
accessory_N=sum((can_grow<(700*.99)).*(can_grow>0));


labels={'Carbon','Nitrogen'};
counts=[core_C,accessory_C; core_N,accessory_N];
C=flip(plasma(4));
X=categorical(labels);
X=reordercats(X,labels);

b=bar(X,counts,'stacked','FaceColor','flat')
for I=1:length(b)
    b(I).CData=0* b(I).CData+C(I,:);
end
legend({'Core ','Accessory'},'Location','northeast')
ylabel('Count')
title('Core and Accessory Sole C/N Sources')
set(gca,'FontName','Arial')
set(gca,'FontSize',14)

%%
clc
clearvars
close all
%another tsne showing no separation by source
load rxn_strain_matrix.mat
load strain_list.mat



x=tsne(1*rxn_strain_matrix','Perplexity',5);

C=flip(plasma(11));

scatter(x(1:344,1),x(1:344,2),30,C(9,:),'filled')
hold on
scatter(x(345:end,1),x(345:end,2),30,C(3,:),'filled')

title('t-SNE of Reaction Presence')
xlabel('Dimension One (Unitless)')
ylabel('Dimension Two (Unitless)')
legend({'Isolates','MAGs'},'Location','northwest')

set(gca,'FontName','Arial')
set(gca,'FontSize',14)