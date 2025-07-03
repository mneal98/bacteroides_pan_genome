%function that gives active exchanges and their bounds.
function [env]=environment(model)
    ex_ind=find(startsWith(model.rxns,'EX_'));
    active_ind=find(model.lb<0);
    active_ex_ind=intersect(ex_ind,active_ind);
    env=[model.rxns(active_ex_ind),num2cell(model.lb(active_ex_ind))];
end