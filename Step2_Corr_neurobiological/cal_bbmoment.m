bbwDir = 'D:\matlab\toolbox\BigBrainWarp';
% load surface
lh = gifti([bbwDir filesep 'template' filesep 'tpl-bigbrain_hemi-L_desc-schaefer2018_7net_400.label.gii']);
rh = gifti([bbwDir filesep 'template' filesep 'tpl-bigbrain_hemi-R_desc-schaefer2018_7net_400.label.gii']);
rh.cdata(rh.cdata>0) = rh.cdata(rh.cdata>0) + 200;
parc_bb = [lh.cdata; rh.cdata];
names = [lh.labels.name rh.labels.name];

% load BigBrain profiles and calculate moments
MP = 65536 - reshape(dlmread([bbwDir '/spaces/tpl-bigbrain/tpl-bigbrain_desc-profiles.txt']),[], 50)';
MPmoments = calculate_moments(MP); % caution will take ~60 minutes to run with 10 parallel workers
writetable(table(MPmoments', parc_bb, 'VariableNames', {'mo', 'yeo'}), ...
            ['D:\Projects\Allo_Scaling\yeo7_moments.csv'])