clc
clear all
%% Problem Settings
lb = [0 0 0 0 0 0 80 100 100 100 100 100];              %lower limit of search space 
ub = [4  4 5 6 9 19 150 250 250 450 450 450];          % upper limit of search space

prob = @Fitness_misfit;    %fitness function

%% Algorithm Parameters
Np = 200;   % number of iteration samples
T = 100 ;     % number of iterations
w=0.73;   %control parameter
c1=1.6;   %control parameter
c2=2.06;   %control parameter

%% Particle Swarmp Optimisation

f = NaN(Np,1);
BestFitIter = NaN(T+1,1);

D = length(lb);

P = repmat(lb,Np,1) + repmat((ub-lb),Np,1).*rand(Np,D);
v = repmat(lb,Np,1) + repmat((ub-lb),Np,1).*rand(Np,D);

for p =1:Np
    f(p) = prob(P(p,:));
end

pbest = P;
f_pbest = f;

[f_gbest,ind] = min(f_pbest);
g_best = P(ind,:);
BestFitIter(1) = f_gbest;


for t = 1:T
    
    for p = 1:Np
        
        v(p,:) = w*v(p,:) + c1*rand(1,D).*(pbest(p,:)-P(p,:)) + c2*rand(1,D).*(g_best - P(p,:));
        
        P(p,:) = P(p,:) + v(p,:);
        
        P(p,:) = max(P(p,:),lb);
        P(p,:) = min(P(p,:),ub);
        
        f(p)=prob(P(p,:));
        
        if f(p) < f_pbest(p)
            
            f_pbest(p) = f(p);
            pbest(p,:) = P(p,:);
            
            if f_pbest(p) < f_gbest
                
                f_gbest = f_pbest(p);
                gbest = pbest(p,:);
                
            end
        end
    end
    
    BestFitIter(t+1)= f_gbest;
    disp(['Iteration' num2str(t) ' : bestfitness = ' num2str(BestFitIter(t+1))]);
end

bestfitness = f_gbest;
bestsol = g_best;


plot(0:T , BestFitIter);
xlabel('Iteration');
ylabel('Best Fitness Value')
title('Particle Swarmp Optimisation');

set(gca,'Fontsize',16,'Fontname','Times New Roman');
set(gcf,'units','centimeters')
% pos = [2, 2, FigWidth, FigHeight]; 
% set(gcf,'Position',pos)
save Result            
