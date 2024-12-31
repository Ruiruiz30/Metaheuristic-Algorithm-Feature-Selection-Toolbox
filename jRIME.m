%[2024]-"Text feature selection using RIME"

% (12/30/2024)
function RIME = jRIME(feat,label,opts)

lb = 0;
ub = 1;
thres = 0.5;
if isfield(opts,'N'), N = opts.N; end
if isfield(opts,'T'), Max_iter = opts.T; end
% Objective function
fun = @jFitnessFunction; 
% Number of dimensions
dim = size(feat,2); 
% Initial 
Rimepop   = zeros(N,dim); 
for i = 1:N
  for d = 1:dim
    Rimepop(i,d) = lb + (ub - lb) * rand();
  end
end  




% initialize position
Best_rime=zeros(1,dim);
Best_rime_rate=inf;%change this to -inf for maximization problems
Lb=lb.*ones(1,dim);% lower boundary 
Ub=ub.*ones(1,dim);% upper boundary
Convergence_curve=zeros(1,Max_iter);
Rime_rates=zeros(1,N);
newRime_rates=zeros(1,N);
W = 5;


for i=1:N
    Rime_rates(1,i)=fun(feat,label,(Rimepop(i,:) > thres),opts);
    
    if Rime_rates(1,i)<Best_rime_rate
        Best_rime_rate=Rime_rates(1,i);
        Best_rime=Rimepop(i,:);
    end
end
Convergence_curve(1) = Best_rime_rate;
newRime_rates = Rime_rates;
it = 2;
% Main loop
while it <= Max_iter
     %Parameters of Eq.(3),(4),(5)
    r1 = (rand-0.5)*2; %[-1,1]
    Sita = (pi*it/(Max_iter/10));
    Beta = (1-round(it*W/Max_iter)/W);
    E =sqrt(it/Max_iter);%Eq.(6)
    RimeFactor = r1 * cos(Sita) * Beta;    
    newRimepop = Rimepop;%Recording new populations
    normalized_rime_rates=normr(Rime_rates);%Parameters of Eq.(7)
    for i=1:N
        for j=1:dim
            
            r1=rand();
            if r1< E
                newRimepop(i,j)=Best_rime(1,j)+RimeFactor*((Ub(j)-Lb(j))*rand+Lb(j));%Eq.(3)
            end
            
            r2=rand();
            if r2<normalized_rime_rates(i)
                newRimepop(i,j)=Best_rime(1,j);%Eq.(7)
            end
        end
    end
    for i=1:N
      
        Flag4ub=newRimepop(i,:)>ub;
        Flag4lb=newRimepop(i,:)<lb;
        newRimepop(i,:)=(newRimepop(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        newRime_rates(1,i)=fun(feat,label,(newRimepop(i,:) > thres),opts);
       
        if newRime_rates(1,i)<Rime_rates(1,i)
            Rime_rates(1,i) = newRime_rates(1,i);
            Rimepop(i,:) = newRimepop(i,:);
            if newRime_rates(1,i)< Best_rime_rate
               Best_rime_rate=Rime_rates(1,i);
               Best_rime=Rimepop(i,:);
            end
        end
    end

    Convergence_curve(it)=Best_rime_rate;
    fprintf('\nIteration %d Best (RIME)= %f',it,Convergence_curve(it))
    it=it+1;
end
% Select features based on selected index
Pos   = 1:dim;
Sf    = Pos((Best_rime > thres) == 1); 
sFeat = feat(:,Sf); 
% Store results
RIME.sf = Sf; 
RIME.ff = sFeat;
RIME.nf = length(Sf);
RIME.c  = Convergence_curve;
RIME.f  = feat;
RIME.l  = label;
end



