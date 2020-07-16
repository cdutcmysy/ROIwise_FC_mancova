# ROIwise_FC_mancova
this function is used for calculate the ROI_wise results from REST
this function is based one spm, so you need to add the spm path in your
% matlab path first.
% this function is used for calculate the ROI_wise results from REST. First, you
% need to move each group in a sub-folder. 
% first you must built a new variates like age/gender as the covariates;
% use example is: [manova_p,sig_p]=ROIwise_FC_mancova(age);
% then select each group's folder, get the mancova_p and sig_p values;
% by YSY & Gallen, Aug, 14, 2018
