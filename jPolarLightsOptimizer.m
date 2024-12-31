% 📜 Polar Lights Optimizer (PLO) Optimization source codes (version 1.0)

%[2024]-"Text feature selection using PLO"

% (12/30/2024)


function PLO=jPolarLightsOptimizer(feat,label,opts)
lb = 0;
ub = 1;
thres = 0.5;
if isfield(opts,'N'), N = opts.N; end
if isfield(opts,'T'), T = opts.T; end
% Objective function
fun = @jFitnessFunction; 
% Number of dimensions
dim = size(feat,2); 

%% Initialization
FEs = 0;
it = 1;
fitness=inf*ones(N,1);
fitness_new=inf*ones(N,1);
Convergence_curve=[];

X=initialization(N,dim,ub,lb);
V=ones(N,dim);
X_new=zeros(N,dim);

for i=1:N
    fitness(i)=fun(feat,label,(X(i,:) > thres),opts);
    FEs=FEs+1;
end

[fitness, SortOrder]=sort(fitness);
X=X(SortOrder,:);
Bestpos=X(1,:);
Bestscore=fitness(1);
Convergence_curve(it)=Bestscore;


%% Main loop
while it <= T
    if it>=T
        break;
    end
    X_sum=sum(X,1);
    X_mean=X_sum/N;
    w1=tansig((FEs/T)^4);
    w2=exp(-(2*FEs/T)^3);
    
    for i=1:N
        
        a=rand()/2+1;
        V(i,:)=1*exp((1-a)/100*FEs);
        LS=V(i,:);

        GS=Levy(dim).*(X_mean-X(i,:)+(lb+rand(1,dim)*(ub-lb))/2);
        X_new(i,:)=X(i,:)+(w1*LS+w2*GS).*rand(1,dim);
    end
    
    E =sqrt(FEs/T);
    A=randperm(N);
    for i=1:N
        for j=1:dim
            if (rand<0.05) && (rand<E)
                X_new(i,j)=X(i,j)+sin(rand*pi)*(X(i,j)-X(A(i),j));
            end
        end
        Flag4ub=X_new(i,:)>ub;
        Flag4lb=X_new(i,:)<lb;
        X_new(i,:)=(X_new(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        fitness_new(i)=fun(feat,label,(X_new(i,:) > thres),opts);
        FEs=FEs+1;
        if fitness_new(i)<fitness(i)
            X(i,:)=X_new(i,:);
            fitness(i)=fitness_new(i);
        end
    end
    [fitness, SortOrder]=sort(fitness);
    X=X(SortOrder,:);
    if fitness(1)<Bestscore
        Bestpos=X(1,:);
        Bestscore=fitness(1);
    end
    if it>=T
        break;
    end
    it = it + 1;
    
    Convergence_curve(it)=Bestscore;
    fprintf('\nIteration %d Best (PLO)= %f',it,Convergence_curve(it))
    Best_pos=Bestpos;
end
best_fun = Bestscore;
 % Select features based on selected index
Pos   = 1:dim;
Sf    = Pos((Best_pos > thres) == 1); 
sFeat = feat(:,Sf); 
% Store results
PLO.sf = Sf; 
PLO.ff = sFeat;
PLO.nf = length(Sf);
PLO.c  = Convergence_curve;
PLO.f  = feat;
PLO.l  = label;
end

function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);
step=u./abs(v).^(1/beta);
o=step;
end