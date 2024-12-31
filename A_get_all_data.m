%---------------------------------------------------------------------------------------------------------------------------%
% Author: Jinrui Zhang

% Last update: 12/31/2024

% E-Mail: ruriuiz@gmail.com

% The useful data from running this code is as follows! 

% A_Fitness_next_sum      ->     Best fitness.                  
% A_Acc_next_sum          ->     Best accuracy.                 
% A_Friedman_Fitness      ->     Friedman rank of fitness.      
% A_Friedman_Acc          ->     Friedman rank of accuracy.     
% A_Fmeasure_next_sum     ->     Best F-measure.                
% A_Precision_next_sum    ->     Best precision.   

% After use of code, please users cite to the main paper: doi.org/10.1007/s10586-024-04901-w

% The toolbox is improve on the works of JingweiToo in github.

%---------------------------------------------------------------------------------------------------------------------------%

clc
clear all
close all

%% Setting parameters
ho = 0.2;         % Ratio of validation data
opts.k = 5;       % Number of k in K-nearest neighbor
opts.N  = 30;     % Number of solutions
opts.T  = 100;    % Maximum number of iterations
Num_Algorithm = 11;        % Number of algorithms.
Num_dataset = 21;     % Number of datasets.
runs = 30;        % Number of run                                         
MAX_runs = 1;     % The maximum number of runs of the entire code

% Opening up array space
Best_Fitness = zeros(Num_dataset,Num_Algorithm);  
Best_Acc = zeros(Num_dataset,Num_Algorithm);   
Mean_Fitness = zeros(Num_dataset,Num_Algorithm);   
Mean_Acc = zeros(Num_dataset,Num_Algorithm);   
Std_Fitness = zeros(Num_dataset,Num_Algorithm);   
Std_Acc = zeros(Num_dataset,Num_Algorithm);     
Mean_time= zeros(Num_dataset,Num_Algorithm); 
Mean_Fmeasure = zeros(Num_dataset,Num_Algorithm); 
Mean_Precision = zeros(Num_dataset,Num_Algorithm); 
% Transition
Fit_sum = zeros(Num_dataset*3,Num_Algorithm); 
Acc_sum = zeros(Num_dataset*3,Num_Algorithm);
Fmeasure_sum = zeros(Num_dataset,Num_Algorithm);
Precision_sum = zeros(Num_dataset,Num_Algorithm);
% Best!
A_Fitness_next_sum = zeros(Num_dataset*3,Num_Algorithm); 
A_Acc_next_sum = zeros(Num_dataset*3,Num_Algorithm);
A_Fmeasure_next_sum = zeros(Num_dataset,Num_Algorithm);
A_Precision_next_sum = zeros(Num_dataset,Num_Algorithm);


Accuracy = zeros(1,Num_Algorithm);    

all_dataF = zeros(Num_Algorithm*Num_dataset,runs);
all_dataA = zeros(Num_Algorithm*Num_dataset,runs);

%% Main
for a = 1:Num_dataset    
    [dataset_name,feat,label] = get_dataset_name(a);
    Best_Fitness = zeros(Num_Algorithm,runs);  
    BestAcc = zeros(Num_Algorithm,runs);
    BestFmeasure = zeros(Num_Algorithm,runs);
    BestPrecision = zeros(Num_Algorithm,runs);
    Time = zeros(Num_Algorithm,runs);
    for p =1:MAX_runs 
        for i = 1:Num_Algorithm  % Run algorithm.
            best_Fitness=zeros(1,runs);  
            best_Acc=zeros(1,runs);      
            best_Fmeasure=zeros(1,runs);
            best_Precision=zeros(1,runs);
            Computational_time=zeros(1,runs); 
            for j = 1:runs  
                disp(" ");
                disp("The dataset for the run are: "+dataset_name+", run "+num2str(p)+" . Run Algorithm "+num2str(i)+" . Running for the "+num2str(j)+" time.") 
                [Algorithm_name,legend_name] = get_Name(i);
                HO = cvpartition(label,'HoldOut',ho);   
                opts.Model = HO;
                FS = jfs(Algorithm_name,feat,label,opts);   
                sf_idx = FS.sf;                                 
                [Acc,Fmeasure,Precision] = Knn_all(feat(:,sf_idx),label,opts);
                best_Fitness(j) = FS.c(length(FS.c));
                best_Acc(j) = Acc;
                best_Fmeasure = Fmeasure;
                best_Precision = Precision;
                Computational_time(j) = FS.t;
            end
            %% Calculate the optimal value, mean, standard deviation
            Best_Fitness(i,:)=sort(best_Fitness);  
            BestAcc(i,:)=sort(best_Acc);
            BestFmeasure(i,:)=sort(best_Fmeasure);
            BestPrecision(i,:)=sort(best_Precision);

            Time(i,:)=sort(Computational_time);

            Best_Fitness(a,i) = Best_Fitness(i,1);   
            Best_Acc(a,i) = BestAcc(i,1); 

            Mean_Fitness(a,i) = mean(Best_Fitness(i,:));  
            Mean_Acc(a,i) = mean(BestAcc(i,:));
            Mean_Fmeasure(a,i) = mean(BestFmeasure(i,:));
            Mean_Precision(a,i) = mean(BestPrecision(i,:));

            Std_Fitness(a,i) = std(Best_Fitness(i,:));   
            Std_Acc(a,i) = std(BestAcc(i,:));

            Mean_time(a,i) = mean(Time(i,:));       
        end
        %% Storing data
        Fit_sum(3*a-2,:)=Best_Fitness(a,:);
        Fit_sum(3*a-1,:) =Mean_Fitness(a,:);
        Fit_sum(3*a,:) = Std_Fitness(a,:);
        Acc_sum(3*a-2,:)=Best_Acc(a,:);
        Acc_sum(3*a-1,:) =Mean_Acc(a,:);
        Acc_sum(3*a,:) = Std_Acc(a,:);
        Fmeasure_sum(a,:) = Mean_Fmeasure(a,:);
        Precision_sum(a,:) = Mean_Precision(a,:);

        %% Get better data
        if p==1 %When p is equal to 1, the data is stored directly
            A_Fitness_next_sum(3*a-2,:) =Fit_sum(3*a-2,:);
            A_Fitness_next_sum(3*a-1,:) =Fit_sum(3*a-1,:);
            A_Fitness_next_sum(3*a,:) =Fit_sum(3*a,:);
            A_Acc_next_sum(3*a-2,:) =Acc_sum(3*a-2,:);
            A_Acc_next_sum(3*a-1,:) =Acc_sum(3*a-1,:);
            A_Acc_next_sum(3*a,:) =Acc_sum(3*a,:);
            A_Fmeasure_next_sum(a,:) = Fmeasure_sum(a,:);
            A_Precision_next_sum(a,:) = Precision_sum(a,:);
            
            %% Used to output all the data of independent runs, where the columns represent the algorithm, and the rows represent the number of times the algorithm was run plus the number of functions
            for j = 1:Num_Algorithm    
                all_dataF((a-1)*Num_Algorithm+j,:)= Best_Fitness(j,:); %all_data used to calculate the Friedman ranking.
                all_dataA((a-1)*Num_Algorithm+j,:)= BestAcc(j,:);
            end
        else
            %% Updating data
            % Fitness
            for i = 1:Num_Algorithm
                if i==1 && Fit_sum(3*a-2,1)<A_Fitness_next_sum(3*a-2,1) && Fit_sum(3*a-1,1)<A_Fitness_next_sum(3*a-1,1) &&  Fit_sum(3*a,1)<A_Fitness_next_sum(3*a,1)
                    A_Fitness_next_sum(3*a-2,i) =Fit_sum(3*a-2,i);
                    A_Fitness_next_sum(3*a-1,i) =Fit_sum(3*a-1,i);
                    A_Fitness_next_sum(3*a,i) =Fit_sum(3*a,i);
                    all_dataF((a-1)*Num_Algorithm+i,:)= Best_Fitness(i,:);
                elseif i~=1
                    if Fit_sum(3*a-1,i)>A_Fitness_next_sum(3*a-1,i) &&  Fit_sum(3*a,i)>A_Fitness_next_sum(3*a,i) &&  Fit_sum(3*a-2,i)>A_Fitness_next_sum(3*a-2,i)
                        A_Fitness_next_sum(3*a-2,i) =Fit_sum(3*a-2,i);
                        A_Fitness_next_sum(3*a-1,i) =Fit_sum(3*a-1,i);
                        A_Fitness_next_sum(3*a,i) =Fit_sum(3*a,i);
                        all_dataF((a-1)*Num_Algorithm+i,:)= Best_Fitness(i,:);
                    end
                end
            end
            % Acc
            for i = 1:Num_Algorithm
                if i==1 && Acc_sum(3*a-2,1)<A_Acc_next_sum(3*a-2,1) && Acc_sum(3*a-1,1)<A_Acc_next_sum(3*a-1,1) &&  Acc_sum(3*a,1)<A_Acc_next_sum(3*a,1)
                    A_Acc_next_sum(3*a-2,i) =Acc_sum(3*a-2,i);
                    A_Acc_next_sum(3*a-1,i) =Acc_sum(3*a-1,i);
                    A_Acc_next_sum(3*a,i) =Acc_sum(3*a,i);
                    all_dataA((a-1)*Num_Algorithm+i,:)= BestAcc(i,:);
                elseif i~=1
                    if Acc_sum(3*a-1,i)>A_Acc_next_sum(3*a-1,i) &&  Acc_sum(3*a,i)>A_Acc_next_sum(3*a,i) &&  Acc_sum(3*a-2,i)>A_Acc_next_sum(3*a-2,i)
                        A_Acc_next_sum(3*a-2,i) =Acc_sum(3*a-2,i);
                        A_Acc_next_sum(3*a-1,i) =Acc_sum(3*a-1,i);
                        A_Acc_next_sum(3*a,i) =Acc_sum(3*a,i);
                        all_dataA((a-1)*Num_Algorithm+i,:)= BestAcc(i,:);
                    end
                end
            end
            % F-measure
            for i = 1:Num_Algorithm
                if i==1 && Fmeasure_sum(a,1)<A_Fmeasure_next_sum(a,1)
                    A_Fmeasure_next_sum(a,i) =Fmeasure_sum(a,i);
                elseif i~=1
                    if Fmeasure_sum(a,i)>A_Fmeasure_next_sum(a,i) 
                        A_Fmeasure_next_sum(a,i) =Fmeasure_sum(a,i);
                    end
                end
            end
            % Precision
            for i = 1:Num_Algorithm
                if i==1 && Precision_sum(a,1)<A_Precision_next_sum(a,1)
                    A_Precision_next_sum(a,i) =Precision_sum(a,i);
                elseif i~=1
                    if Precision_sum(a,i)>A_Precision_next_sum(a,i) 
                        A_Precision_next_sum(a,i) =Precision_sum(a,i);
                    end
                end
            end
        end
    end
end

all_data_Fitness = all_dataF';
all_data_acc = all_dataA';
%% Calculating the fitness Friedman
RawData = all_data_Fitness;
% Finding Rank of the problems
for j = 1:Num_dataset
    for i = 1:runs
        RankOfTheProblems(i,:) = tiedrank(RawData(i,1+(j-1)*Num_Algorithm:j*Num_Algorithm));
    end
    % Taking average of the rank of the problems
    AvgOfRankOfProblems = mean(RankOfTheProblems);
    SquareOfTheAvgs = AvgOfRankOfProblems .* AvgOfRankOfProblems;
    SumOfTheSquares = sum(SquareOfTheAvgs);
    FfStats = (12*runs/(Num_Algorithm*(Num_Algorithm+1))) * (SumOfTheSquares - ((Num_Algorithm*(Num_Algorithm+1)^2)/4));
    A_Friedman_Fitness(j,:) = AvgOfRankOfProblems;
end
Friedman_Mean_Fitness = mean(A_Friedman_Fitness);
[a,b]=sort(Friedman_Mean_Fitness);
for i = 1:size(b,2)
    Friedman_FinalM_Fitness(b(i)) = i;
end
disp(" ")
disp("%%%%%%%%%%%%%%%%%%%%%%%Fitness Ranking - The value of the output is%%%%%%%%%%%%%%%%%%%%%%");
disp("Ranking of individual algorithms on different functions：");
disp(A_Friedman_Fitness)
disp("Average rank：");
disp(Friedman_Mean_Fitness);
disp("Final average rank：");
disp(Friedman_FinalM_Fitness);
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");

%% Calculating accuracy Friedman
RawData = all_data_acc.*(-1);
% Finding Rank of the problems
for j = 1:Num_dataset
    for i = 1:runs
        RankOfTheProblems(i,:) = tiedrank(RawData(i,1+(j-1)*Num_Algorithm:j*Num_Algorithm));
    end
    % Taking average of the rank of the problems
    AvgOfRankOfProblems = mean(RankOfTheProblems);
    SquareOfTheAvgs = AvgOfRankOfProblems .* AvgOfRankOfProblems;
    SumOfTheSquares = sum(SquareOfTheAvgs);
    FfStats = (12*runs/(Num_Algorithm*(Num_Algorithm+1))) * (SumOfTheSquares - ((Num_Algorithm*(Num_Algorithm+1)^2)/4));
    % Display the results
    %formatSpec = 'Friedman statistic is %4.2f and \n ';
    %fprintf(formatSpec,FfStats);
    %disp('Average of the ranks obtained in all problems');
    %disp(AvgOfRankOfProblems)
    A_Friedman_Acc(j,:) = AvgOfRankOfProblems;
    %y(j,:)= tiedrank(AvgOfRankOfProblems);
end
Friedman_Mean_Acc = mean(A_Friedman_Acc);
[a,b]=sort(Friedman_Mean_Acc);
for i = 1:size(b,2)
    Friedman_FinalM_Acc(b(i)) = i;
end
disp("%%%%%%%%%%%%%%%%%%%%%%%Accuracy Ranking - The value of the output is%%%%%%%%%%%%%%%%%%%%%%%");
disp("Ranking of individual algorithms on different functions：");
disp(A_Friedman_Acc)
disp("Average rank：");
disp(Friedman_Mean_Acc);
disp("Final average rank：");
disp(Friedman_FinalM_Acc);
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");





