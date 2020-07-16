% this function is based one spm, so you need to add the spm path in your
% matlab path first.
% this function is used for calculate the roi-wise FC analysis results from REST. First, you
% need to move each group FC matrix in a sub-folder.
% use example is: [manova_p]=ROIwise_FC_manova_nets;
% then select each group's folder,
% For YSY, by Gallen, Aug, 14, 2018
function [manova_p,manova_R2,manova_F,delta_value]=ROIwise_FC_mancova_nets(covariates)
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
    g3_data(:,:,i)=load([path_g1 filesep file_g1(i).name]);
end

for j = 1:length(file_g2)
    g4_data(:,:,j)=load([path_g2 filesep file_g2(j).name]);
end
num1 = size(g1_data,3);num2 = size(g3_data,3);
for l = 1:num1
    treatment{l,:} = treatment1{1};
end
for l = (num1+1):(num1+num2)
    treatment{l,:} = treatment2{1};
end
prompt = {'Enter a value of each nets number'};
titles = 'Each nets number';
definput = {'4 5 4 3 2'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,titles,[1 40],definput,opts);
netnumber = str2num(answer{1});
if size(g1_data,1)~=sum(netnumber)
    error('The nets sum number is not equal to the input matrix, please check it!')
end
% obtain the networks pair and compare times
comp_time = length(netnumber)+length(netnumber)*((length(netnumber)-1))/2;
aa = ones(length(netnumber));
a=tril(aa);
[s(:,1),s(:,2)]=find(a==1);
paired_nets = s;
% get the average FC in each paried networks
mm = 1;
for i = 1:length(netnumber)
    mm(i+1) = mm(i) + netnumber(i);
    space{i} = mm(i):(mm(i+1)-1);
end
for i = 1:comp_time
    m = s(i,1);n = s(i,2);
    for j = 1:size(g1_data,3)
        tmp = g1_data(space{m},space{n},j);
        tmp(tmp>10)=0;
        g1_net_data(j,i) = mean(tmp(tmp~=0));
    end
    for j = 1:size(g2_data,3)
        tmp = g2_data(space{m},space{n},j);
        tmp(tmp>10)=0;
        g2_net_data(j,i) = mean(tmp(tmp~=0));
    end
    for j = 1:size(g3_data,3)
        tmp = g3_data(space{m},space{n},j);
        tmp(tmp>10)=0;
        g3_net_data(j,i) = mean(tmp(tmp~=0));
    end
    for j = 1:size(g4_data,3)
        tmp = g4_data(space{m},space{n},j);
        tmp(tmp>10)=0;
        g4_net_data(j,i) = mean(tmp(tmp~=0));
    end
end
data1 = cat(1,g1_net_data,g3_net_data);
data2 = cat(1,g2_net_data,g4_net_data);
data(:,:,1) = data1;data(:,:,2)=data2;
delta_value = data(:,:,1)-data(:,:,2);
% paired t-test for each paired nets
for n = 1:size(data,2)
    Y = squeeze(data(:,n,:));
    y = Y;
    covariates = covariates;
    [b,bint,r,rint,stats] = regress(y(:,1),[ones(length(y),1) covariates]);
    Y_new = y(:,1)-covariates*b(2:end);
    [b,bint,r,rint,stats] = regress(y(:,2),[ones(length(y),1) covariates]);
    Y_new(:,2) = y(:,2)-covariates*b(2:end);
    
    t = table(treatment,Y_new(:,1),Y_new(:,2),'VariableNames',{'treatment','pre','post'});
    Meas = table([1 2]','VariableNames',{'timepoint'});
    rm = fitrm(t,'pre-post~treatment','WithinDesign',Meas);
    [manovatbl,A,C,D] = manova(rm);
    manova_p(n,1) = manovatbl.pValue(5);
    manova_R2(n,1) = manovatbl.RSquare(5);
    manova_F(n,1) = manovatbl.F(5);
end
[~,FDR_p] = FDR(manova_p,0.05); % FDR correction
save ROIwise_mancova_nets_age_d_noGS.mat
