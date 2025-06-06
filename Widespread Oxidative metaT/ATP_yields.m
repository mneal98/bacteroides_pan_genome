clc
clearvars
close all

%here we examine ATP yields to see if all strains have effective oxygen and
%nitrite utilization pathways. 
load pan_model.mat
load rxn_strain_matrix.mat
rxn_strain_matrix((strcmp(model.rxns,'NO2tex')),:)=1;

ATP_yield=zeros(700,3); %anaerobic, with o2, with no2
for I=1:700
    strain_model=removeRxns(model,model.rxns(~rxn_strain_matrix(:,I)));
    strain_model=changeObjective(strain_model,'ATPM');

    FBA=optimizeCbModel(strain_model);
    ATP_yield(I,1)=FBA.f;

    strain_model.lb(strcmp(strain_model.rxns,'EX_o2_e'))=-1000;
    FBA=optimizeCbModel(strain_model);
    ATP_yield(I,2)=FBA.f;
    strain_model.lb(strcmp(strain_model.rxns,'EX_o2_e'))=0;

    strain_model.lb(strcmp(strain_model.rxns,'EX_no2_e'))=-1000;
    FBA=optimizeCbModel(strain_model);
    ATP_yield(I,3)=FBA.f;

    if mod(I,50)==0
        fprintf(num2str(I))
        fprintf('\n')
    end
end
ATP_yield=ATP_yield/10;
improvement=[ATP_yield(:,2)./ATP_yield(:,1),ATP_yield(:,3)./ATP_yield(:,1)];