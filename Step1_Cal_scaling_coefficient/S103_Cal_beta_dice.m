clear
%% dice
cd D:\Projects\Allo_Scaling\
load Data\Beta\beta_deg_hcp.mat
FCN_sig_HCP = find(beta_all.FC.valid_beta~=0);
MCN_sig_HCP = find(beta_all.GMV.valid_beta~=0);
sum(beta_all.SC.p_beta<=beta_all.SC.sig_thr)
sum(beta_all.FC.p_beta<=beta_all.FC.sig_thr)
sum(beta_all.GMV.p_beta<=beta_all.GMV.sig_thr)

load Data\Beta\beta_deg_slem.mat
FCN_sig_SLEM = find(beta_all.FC.valid_beta~=0);
MCN_sig_SLEM = find(beta_all.GMV.valid_beta~=0);
sum(beta_all.SC.p_beta<=beta_all.SC.sig_thr)
sum(beta_all.FC.p_beta<=beta_all.FC.sig_thr)
sum(beta_all.GMV.p_beta<=beta_all.GMV.sig_thr)


gretna_dice_coefficient(FCN_sig_HCP,FCN_sig_SLEM)
gretna_dice_coefficient(MCN_sig_HCP,MCN_sig_SLEM)