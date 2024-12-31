%---------------------------------------------------------------------------------------------------------------------------%
% Author: Jinrui Zhang

% Last update: 12/31/2024

% E-Mail: ruriuiz@gmail.com

% The useful data from running this code is as follows! 

% After use of code, please users cite to the main paper: doi.org/10.1007/s10586-024-04901-w

% The toolbox is improve on the works of JingweiToo in github.

%---Input-------------------------------------------------------------
% feat   : Feature vector matrix (Instances x Features)
% label  : Label matrix (Instances x 1)
% opts   : Parameter settings
% opts.N : Number of solutions / population size (* for all methods)
% opts.T : Maximum number of iterations (* for all methods)
% opts.k : Number of k in k-nearest neighbor

%---Output------------------------------------------------------------
% FS    : Feature selection model (It contains several results)
% FS.sf : Index of selected features
% FS.ff : Selected features
% FS.nf : Number of selected features
% FS.c  : Convergence curve
% Acc   : Accuracy of validation model

clc
clear all;
close all;

% Setting parameters
ho = 0.2;           % Ratio of validation data
opts.k = 5;         % Number of k in K-nearest neighbor
opts.N  = 10;       % number of solutions
opts.T  = 100;      % maximum number of iterations
Num_Algorithm = 1;  % Number of algorithms.
Num_dataset = 21;   % Number of datasets.
[dataset,feat,label] = get_Dataset(Num_dataset);
Accuracy = zeros(1,Num_Algorithm);                      
Precision = zeros(1,Num_Algorithm);                      
F = zeros(1,Num_Algorithm);                                 
Computational_time = zeros(1,Num_Algorithm);         
size = zeros(1,Num_Algorithm);
c = categorical({});                        
color = ["r","c","k","y","m","g","b","[1 0.5 0]","[0 0.4470 0.7410]","[0.3010 0.7450 0.9330]"];
for j=1:Num_Algorithm
    [Algorithm_name,legend_name] = get_Algorithm(j);
    c(j) = legend_name;
    HO = cvpartition(label,'HoldOut',ho);         
    opts.Model = HO;
    FS     = jfs(Algorithm_name,feat,label,opts);   
    sf_idx = FS.sf;                                
    [Acc,Fmeasure,Precision] = Knn_all(feat(:,sf_idx),label,opts);
    Accuracy(1,j) = Acc;
    Computational_time(1,j) = FS.t;
    size(1,j) = FS.nf;
    if j==1
        plot(FS.c,'Color',color(j),'LineWidth',3);
        hold on;
    else
        plot(FS.c,'Color',color(j),'LineWidth',1.5);
        hold on;
    end

end
%Fitness
grid on;
legend(c);
set(legend,'NumColumns',3);
xlabel('Number of Iterations');
ylabel('Fitness Value');
title(dataset);
% Accuracy
figure(2);
color2 = ["[1 0.3098 0.3098]","[0.1490 0.5176 0.9412]","[0.51765 0.43922 1]","[0.6863 0.9804 0.5765]","[1 1 0]","[1 0.7098 0.9608]","[0.3922 1 0.2392]","[0.9882 0.5333 0.3647]","[0.8196 0.6902 0.1725]","[0.4627 0.7176 0.6980]"];
for k = 1:Num_Algorithm
    bar(c(k),Accuracy(k),'FaceColor',color2(k));
    hold on;
end
ylim([0.6,1.0]);
grid on;
ylabel('Accuracy');
title(dataset);
% Feature size
figure(3);
for k = 1:Num_Algorithm
    bar(c(k),size(k),'FaceColor',color2(k));
    hold on;
end
grid on;
ylabel('Feature Size');
title(dataset);

