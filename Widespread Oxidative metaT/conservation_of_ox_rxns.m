%Here are the main oxidative pathways. We show that these are present in 
%the vast majority of strains.

clc
clearvars
close all

load rxn_strain_matrix.mat
load pan_model.mat


% cytochrome oxidase, nitrite reductase,   nadh dehydrogenase, ATP synthase
ox_rxns={'CYTBD2pp',   'NTRIR4pp',        'NQR',              'ATPS',...
         'SPODMpp',      'CAT'            ,'GND',...
         'CS','ACONTa','ICDHyr','AKGDH','SUCOAS','FRD2rpp','FUM','MDH'};
            %superoxide dismutase, catalase, oxidative pentose phosphate, then TCA
for I=1:length(ox_rxns)
    percent_strains(I,1)=sum(rxn_strain_matrix(strcmp(model.rxns,ox_rxns{I}),:));
end
percent_strains=percent_strains/700*100;
percent_strains=round(percent_strains,1)
%and thus these rxns are present in 97+% of strains. 