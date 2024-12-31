%----------------------------------------------------------
% Crayfish Optimization Algorithm
% https://doi.org/10.1007/s10462-023-10567-4
%
%[2024]-"Text feature selection using PLO"
%
% (12/30/2024)
%----------------------------------------------------------
function COA = jCrayfishOptimizationAlgorithm(feat,label,opts)
lb = 0;
ub = 1;
thres = 0.5;
if isfield(opts,'N'), N = opts.N; end
if isfield(opts,'T'), T = opts.T; end
% Objective function
fun = @jFitnessFunction; 
% Number of dimensions
dim = size(feat,2); 
% Initial 
X   = zeros(N,dim); 
for i = 1:N
  for d = 1:dim
    X(i,d) = lb + (ub - lb) * rand();
  end
end  
%% Define Parameters
cuve_f=zeros(1,T); 
global_Cov = zeros(1,T);
Best_fitness = inf;
best_position = zeros(1,dim);
fitness_f = zeros(1,N);
t=1; 
for i=1:N
   fitness_f(i) =  fun(feat,label,(X(i,:) > thres),opts); %Calculate the fitness value of the function
   if fitness_f(i)<Best_fitness
       Best_fitness = fitness_f(i);
       best_position = X(i,:);
   end
      cuve_f(t)=Best_fitness;
end
global_position = best_position; 
global_fitness = Best_fitness;


p=[1:360];
while(t<=T)
    C = 2-(t/T); %Eq.(7)
    temp = rand*15+20; %Eq.(3)
    xf = (best_position+global_position)/2; % Eq.(5)
    Xfood = best_position;
    
    
    for i = 1:N
        if temp>30                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
           
            %% summer resort stage
            if rand<0.5
                Xnew(i,:) = X(i,:)+C*rand(1,dim).*(xf-X(i,:)); %Eq.(6)
            else
            %% competition stage
                for j = 1:dim
                    z = round(rand*(N-1))+1;  %Eq.(9)
                    Xnew(i,j) = X(i,j)-X(z,j)+xf(j);  %Eq.(8)
                end
            end
            
        else
            %% foraging stage
            F1 = fun(feat,label,(Xfood > thres),opts);
            P = 3*rand*fitness_f(i)/F1; %Eq.(4)
            if P>2   % The food is too big
                 Xfood = exp(-1/P).*Xfood;   %Eq.(12)
                for j = 1:dim
                    Xnew(i,j) = X(i,j)+cos(2*pi*rand)*Xfood(j)*p_obj(temp,t,T)-sin(2*pi*rand)*Xfood(j)*p_obj(temp,t,T); %Eq.(13)
                end
            else
                Xnew(i,:) = (X(i,:)-Xfood)*p_obj(temp,t,T)+p_obj(temp,t,T).*rand(1,dim).*X(i,:); %Eq.(14)
            end
            
        end
    end
    %% boundary conditions
    for i=1:N
        for j =1:dim
            if length(ub)==1
                Xnew(i,j) = min(ub,Xnew(i,j));
                Xnew(i,j) = max(lb,Xnew(i,j));
            else
                Xnew(i,j) = min(ub(j),Xnew(i,j));
                Xnew(i,j) = max(lb(j),Xnew(i,j));
            end
        end
    end
   
    global_position = Xnew(1,:);
    global_fitness = fun(feat,label,(global_position > thres),opts);
    for i =1:N
         %% Obtain the optimal solution for the updated population
        new_fitness = fun(feat,label,(Xnew(i,:) > thres),opts);
        if new_fitness<global_fitness
                 global_fitness = new_fitness;
                 global_position = Xnew(i,:);
                
        end
        %% Update the population to a new location
        if new_fitness<fitness_f(i)
             fitness_f(i) = new_fitness;
             X(i,:) = Xnew(i,:);
             if fitness_f(i)<Best_fitness
                 Best_fitness=fitness_f(i);
                 best_position = X(i,:);
             end
        end
        cuve_f(t) = Best_fitness;

    end
    global_Cov(t) = global_fitness;
    t = t+1;

end
 best_fun = Best_fitness;
 % Select features based on selected index
Pos   = 1:dim;
Sf    = Pos((best_position > thres) == 1); 
sFeat = feat(:,Sf); 
% Store results
COA.sf = Sf; 
COA.ff = sFeat;
COA.nf = length(Sf);
COA.c  = cuve_f;
COA.f  = feat;
COA.l  = label;
end
function y = p_obj(x,t,T)   %Eq.(4)
C=0.2;
Q=3;

y = C*(1/(sqrt(2*pi)*Q))*exp(-(x-25).^2/(2*Q.^2));
end