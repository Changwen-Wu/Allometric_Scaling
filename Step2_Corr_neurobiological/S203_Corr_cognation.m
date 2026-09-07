cd D:\Projects\Allo_Scaling

load Data\Beta\beta_deg_hcp.mat
load Data\Spin_rotation\Schaefer400_rotation.mat

Permutation_times = 10000;

neurosynth_tbl = readtable("Data\Neurosynth\atl-schaefer2018_res-400_neurosynth.csv");
activation_score = table2array(neurosynth_tbl(:,2:end));
terms = neurosynth_tbl.Properties.VariableNames(2:end);

% mode = {'SC', 'FC', 'GMV'};
mode = {'FC'};
P = zeros(length(mode),size(activation_score,2));
R = zeros(length(mode),size(activation_score,2));

thr = zeros(1,length(mode));

for imod = 1 : length(mode)
    beta = beta_all.(mode{imod}).beta;

    for iterm = 1 : size(activation_score,2)
        if mod(iterm,20) == 0; disp(num2str(iterm/size(activation_score,2))); end
        data = activation_score(:,iterm);
        if size(data, 1) == 1; data = data.'; end
        if size(beta, 1) == 1; beta = beta.'; end

        [P(imod,iterm), ~] = perm_sphere_p_wcw(beta, data, perm_id, 'spearman');
        [R(imod,iterm),~] = corr(beta,data,'type','Spearman','rows','pairwise');
    end
    thr(imod) = gretna_FDR(P(imod,:),0.05);
end

save Result\S203_Corr_cognation\corr_cognation_hcp.mat P R thr terms

fc_p = P; fc_r = R; fc_thr = thr;
sig_pos_term = terms(and(fc_p<=fc_thr, fc_r>0)); sig_pos_r = fc_r(and(fc_p<=fc_thr, fc_r>0));
sig_neg_term = terms(and(fc_p<=fc_thr, fc_r<0)); sig_neg_r = fc_r(and(fc_p<=fc_thr, fc_r<0));
[sig_pos_r,seq] = sort(sig_pos_r); sig_pos_term = sig_pos_term(seq);
[sig_neg_r,seq] = sort(sig_neg_r); sig_neg_term = sig_neg_term(seq);