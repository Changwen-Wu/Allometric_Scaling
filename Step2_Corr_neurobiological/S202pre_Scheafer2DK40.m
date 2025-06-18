clear; clc;
cd D:\Projects\Allo_Scaling

%% left hemisphere
% map1
[~, map1, colortable] = read_annotation('lh.Schaefer2018_400Parcels_7Networks_order.annot');
label_map1 = colortable.table(2:end,5);
map1_size = length(unique(label_map1));

% map2
[~, map2, colortable_map2] = read_annotation('lh.DK40_70_fsaverage_164k.annot');
label_map2 = colortable_map2.table(2:end,5);
map2_size = length(unique(label_map2));

Dice_LH = zeros(map1_size, map2_size);

for i = 1 : map1_size
    val1 = label_map1(i);
    for j = 1 : map2_size
        val2 = label_map2(j);
        Dice_LH(i,j) = sum(map1(map2 == val2) == val1)/ sum(map2 == val2);
    end
end

%% right hemisphere
% map1
[~, map1, colortable] = read_annotation('rh.Schaefer2018_400Parcels_7Networks_order.annot');
label_map1 = colortable.table(2:end,5);
map1_size = length(unique(label_map1));

% map2
[~, map2, colortable_map2] = read_annotation('rh.DK40_70_fsaverage_164k.annot');
label_map2 = colortable_map2.table(2:end,5);
map2_size = length(unique(label_map2));

Dice_RH = zeros(map1_size, map2_size);

for i = 1 : map1_size
    val1 = label_map1(i);
    for j = 1 : map2_size
        val2 = label_map2(j);
        Dice_RH(i,j) = sum(map1(map2 == val2) == val1)/ sum(map2 == val2);
    end
end

save Data\Disorder\DK_Sch_overlap.mat Dice_LH Dice_RH
