% this function is based one spm, so you need to add the spm path in your
% matlab path first.
% this function is used for calculate the ROI_wise results from REST. First, you
% need to move each group in a sub-folder. 
% first you must built a new variates like age/gender as the covariates;
% use example is: [manova_p,sig_p]=ROIwise_FC_mancova(age);
% then select each group's folder, get the mancova_p and sig_p values;
% by YSY & Gallen, Aug, 14, 2018
function [manova_p,manova_R2,manova_F,sig_p]=ROIwise_FC_mancova(covariates)
%prepare group1
covariates = covariates;
% covariates = cov;
path_g1 = spm_select(1,'dir','Please selcet group1 pre dir');
path_g2 = spm_select(1,'dir','Please selcet group1 post dir');
prompt = {'Enter the treatment 1 name'};
titles = 'treatment1 name';
definput = {'real'};
opts.Interpreter = 'tex';
treatment1 = inputdlg(prompt,titles,[1 40],definput,opts);

file_g1 = dir([path_g1 filesep '*.txt']);
file_g2 = dir([path_g2 filesep '*.txt']);

for i = 1:length(file_g1)
    g1_data(:,:,i)=load([path_g1 filesep file_g1(i).name]);
end

for j = 1:length(file_g2)
    g2_data(:,:,j)=load([path_g2 filesep file_g2(j).name]);
end

data1(:,:,:,1) = g1_data;
data1(:,:,:,2) = g2_data;

clear g1_data g2_data file_g1 file_g2
% treatment 2 
path_g1 = spm_select(1,'dir','Please selcet group2 pre dir');
path_g2 = spm_select(1,'dir','Please selcet group2 post dir');
prompt = {'Enter the treatment 2 name'};
titles = 'treatment1 name';
definput = {'sham'};
opts.Interpreter = 'tex';
treatment2 = inputdlg(prompt,titles,[1 40],definput,opts);


file_g1 = dir([path_g1 filesep '*.txt']);
file_g2 = dir([path_g2 filesep '*.txt']);

for i = 1:length(file_g1)
    g1_data(:,:,i)=load([path_g1 filesep file_g1(i).name]);
end

for j = 1:length(file_g2)
    g2_data(:,:,j)=load([path_g2 filesep file_g2(j).name]);
end

data2(:,:,:,1) = g1_data;
data2(:,:,:,2) = g2_data;

num1 = size(data1,3);num2 = size(data2,3);
for l = 1:num1
treatment{l,:} = treatment1{1};
end
for l = (num1+1):(num1+num2)
treatment{l,:} = treatment2{1};
end
data = cat(3,data1,data2);

for m = 1:size(data,1)
    for n = 1:size(data,2)
        if m~=n
        Y = squeeze(data(m,n,:,:));
        y = Y;
        y(~isfinite(y))=0;
        %     covariates = aa';
        [b,bint,r,rint,stats] = regress(y(:,1),[ones(length(y),1) covariates]);
        Y_new = y(:,1)-covariates*b(2:end);
        [b,bint,r,rint,stats] = regress(y(:,2),[ones(length(y),1) covariates]);
        Y_new(:,2) = y(:,2)-covariates*b(2:end);
        
        t = table(treatment,Y_new(:,1),Y_new(:,2),'VariableNames',{'treatment','pre','post'});
        Meas = table([1 2]','VariableNames',{'timepoint'});
        rm = fitrm(t,'pre-post~treatment','WithinDesign',Meas);
        [manovatbl,A,C,D] = manova(rm);
        manova_p(m,n,1) = manovatbl.pValue(5);
        manova_R2(m,n,1) = manovatbl.RSquare(5);
        manova_F(m,n,1) = manovatbl.F(5);
        end
    end
end
line_p = tril(manova_p,-1); line_p(line_p==0)=[];
for ll = 1:size(data,3)
   tmp = data(:,:,ll,1);
   data1_line = tril(tmp,-1);data1_line(data1_line==0)=[];
   pre_data(ll,:) = data1_line;
   tmp = data(:,:,ll,2);
   data1_line = tril(tmp,-1);data1_line(data1_line==0)=[];
   post_data(ll,:) = data1_line;
end
sig_data_pre = pre_data(:,find(line_p<0.05));
sig_data_post = post_data(:,find(line_p<0.05));
sig_data_change = sig_data_post-sig_data_pre;
[~,FDR_p] = FDR(line_p,0.05);
% [~,FDR_p] = FDR(manova_p,0.05);
[sig_p(:,1),sig_p(:,2)] = find(manova_p<0.05);
p_005(:,1) = find(line_p<0.05);
p_005(:,2) = line_p(find(line_p<0.05));
% p_005(:,3) = ROIsequence(find(diff_p<0.05));
% p_005(:,4) = ROIsequence(find(diff_p<0.05),2);
sig_p = p_005;
save ROI_FC_mancova.mat

