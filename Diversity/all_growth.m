clc
clearvars
pan_model=importdata('pan_model.mat');

load rxn_strain_matrix.mat
load strain_list.mat


exchanges=pan_model.rxns(startsWith(pan_model.rxns,'EX_'));
growth_rate=zeros(length(strain_list),length(exchanges));
%REPLACE GLC WITH NH4 for C SOURCE TESTING
% no_C={'pi_e','glc__D_e','nh4_e','HEME','k_e','cobalt2_e','mg2_e','mn2_e','ca2_e','h_e','h2o_e'};
no_C=environment(pan_model);
no_C=no_C(:,1);
no_C=setdiff(no_C,'EX_glc__D_e');
pan_model.lb(startsWith(pan_model.rxns,'EX_'))=0;
%set minerals and such on
rxn_ind=[];
for I=1:length(no_C)
    rxn_ind(I)=find(startsWith(pan_model.rxns,[no_C{I}]));
end
pan_model.lb(rxn_ind)=-10;
%first we find where the pan_model can grow so we can skip anything else in
%the exhcnage list. if it cant do it no subset of rxns could
for I=1:length(exchanges)
    t=pan_model.lb(strcmp(pan_model.rxns,exchanges{I}));
    pan_model.lb(strcmp(pan_model.rxns,exchanges{I}))=-10;
    FBA=optimizeCbModel(pan_model);
    good_ex(I)=FBA.f>.001;
    pan_model.lb(strcmp(pan_model.rxns,exchanges{I}))=t;
end

sum(good_ex)

environment = getEnvironment();
for strain_ind=1:length(strain_list)
restoreEnvironment(environment)
    model=removeRxns(pan_model,pan_model.rxns(~rxn_strain_matrix(:,strain_ind)));
    V=zeros(size(exchanges));

    for met_ind=1:length(exchanges)
        
        if ~good_ex(met_ind)
            continue
        end
        
        t=model.lb(strcmp(model.rxns,exchanges{met_ind}));
        model.lb(strcmp(model.rxns,exchanges{met_ind}))=-10;



        FBA=optimizeCbModel(model,'max');
        V(met_ind)=FBA.f;
        model.lb(strcmp(model.rxns,exchanges{met_ind}))=t;
    end
    
    strain_ind
    growth_rate(strain_ind,:)=V;
end
%remove NaNs
growth_rate(isnan(growth_rate))=0;