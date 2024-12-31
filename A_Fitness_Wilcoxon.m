clear all;
close all;
clc
%% Setting parameters
ho = 0.2;         % Ratio of validation data
opts.k = 5;       % Number of k in K-nearest neighbor
opts.N  = 10;     % number of solutions
opts.T  = 100;    % Maximum number of iterations
Num_Algorithm = 10     ;  % Define how many algorithms need to be run
dataset = 16;     % Number of datasets
runs=30;          % Number of Wilcoxon runs                                         
dat = zeros(Num_Algorithm, runs);
Wilcoxon=zeros(dataset,Num_Algorithm-1);
BestF = zeros(Num_Algorithm,dataset,runs);
for i = 1:dataset
     [dataset_name,feat,label] = get_dataset_name(i);
    for run=1:runs
        for a=1:Num_Algorithm
            [Algorithm_name,legend_name] = get_Name(a);
                HO = cvpartition(label,'HoldOut',ho);      
                opts.Model = HO;
                FS     = jfs(Algorithm_name,feat,label,opts);   %Perform feature selection
                sf_idx = FS.sf;                                 
                [Acc] = jknn(feat(:,sf_idx),label,opts);
            BestF(a,i,run) = FS.c(length(FS.c));
        end
        disp(run)
    end
end

%% Wilcoxon signed-rank test, The first algorithm is compared with other algorithms
for i=1:dataset
    for a=1:Num_Algorithm
        dat(a,:) = BestF(a,i,:);
    end

    Wilcoxon(i,:)=my_pv(['F',num2str(i)], dat, Algorithm_name);
end
