%This shows an example of reconstructing and optimizing a single
%strain-specific model. This should take <1 second.
clc
clearvars

%load in the pan-genome model
pan_model=importdata('pan_model.mat');

%load in the matrix encoding strain-reaction associations
rxn_strain_matrix=importdata('rxn_strain_matrix.mat');

%lets make two strain specific models and compare their growth rates.
strain_141_model=removeRxns( pan_model, pan_model.rxns(rxn_strain_matrix(:,141)==0) );
FBA=optimizeCbModel(strain_141_model);
['Growth rate of strain 141: ',num2str(FBA.f)]

strain_80_model=removeRxns( pan_model, pan_model.rxns(rxn_strain_matrix(:,80)==0) );
FBA=optimizeCbModel(strain_80_model);
['Growth rate of strain 80: ',num2str(FBA.f)]