clear;
cd E:\Projects\Allo_Scaling
%%
load Data\Beta\beta_deg_slem.mat
load Data\Spin_rotation\Schaefer400_rotation.mat
load Data\Neuro_signatures\Neuro_signatures.mat
mode = {'FC'};

corr_signatures = {};

%% hierarchy 1-7
% ind = 1:7;
% label = labels(ind);
% [result] = spin_spearman(signatures,beta_all,ind,perm_id);
% result.label = label; corr_signatures.Hierarchy = result; clear result;


%% neurochemo 8-25

ind = 8:25;
label = labels(ind);
[result] = spin_spearman(signatures,beta_all,ind,perm_id);
result.label = label; corr_signatures.Neurochemo = result; clear result;


%% Gene expression 26-32

ind = 26:46;
label = labels(ind);
[result] = spin_spearman(signatures,beta_all,ind,perm_id);
result.label = label; corr_signatures.Gene = result; clear result;


%% microstructure 47-51

% ind = 47:51;
% label = labels(ind);
% [result] = spin_spearman(signatures,beta_all,ind,perm_id);
% result.label = label; corr_signatures.Microstructure = result; clear result;


% %% Neuropeptide 52-89
% 
% ind = 52:89;
% label = labels(ind);
% [result] = spin_spearman(signatures,beta_all,ind,perm_id);
% result.label = label; corr_signatures.Neuropeptide = result; clear result;


% %% Metabolism 90-94
% 
% ind = 90:94;
% label = labels(ind);
% [result] = spin_spearman(signatures,beta_all,ind,perm_id);
% result.label = label; corr_signatures.Metabolism = result; clear result;


% %% Mitochondrial_features
% load Data\Neuro_signatures\Mitochondrial_features\Mito_feature_sch400.mat data labels
% 
% [result] = spin_spearman(data,beta_all,1:6,perm_id);
% result.label = labels; corr_signatures.Mitcochondrial = result; clear result;

save Result\S201_Corr_signatures\slem_corr_signature.mat corr_signatures mode


function [result] = spin_spearman(signatures,beta_all,ind,perm_id)

if size(signatures, 2) == 400; signatures = signatures.'; end
signatures = signatures(:,ind);

mode = {'FC'};
result = {};

P = zeros(1,size(signatures, 2));
R = zeros(1,size(signatures, 2));

thr = zeros(1,length(mode));

for imod = 1 : length(mode)
    beta = beta_all.(mode{imod}).beta;
    for ifield = 1 : size(signatures, 2)
        data = signatures(:,ifield);

        if size(beta, 1) == 1; beta = beta.'; end

        [P(imod,ifield), ~] = perm_sphere_p_wcw(beta, data, perm_id, 'spearman');
        [R(imod,ifield),~] = corr(beta,data,'type','Spearman','rows','pairwise');
    end
   thr(imod) = gretna_FDR(P(imod,:),0.05);
end

result.P = P; result.R = R; result.thr = thr;

end
