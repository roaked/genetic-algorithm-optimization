% function [AvgCost, BestCost, GlobalBest, BestIt, it] = PSO(model)
% clc
% clear all
% close all
% %% Problem Definition
% model = CreateModel(2);

% model = CreateModel(2);
CostFunction = @(x) PSO_BinPackingCost(x, model);  % Objective Function

nVar = 2*model.n-1;     % Number of Decision Variables
VarSize = [1 nVar];     % Decision Variables Matrix Size

VarMin = 0;     % Lower Bound of Decision Variables
VarMax = 1;     % Upper Bound of Decision Variables


%% PSO Parameters

%MaxIt = 7500

nPop = PSOdata.nPop;           % Population Size (Swarm Size)

w= 0;%PSOdata.w;%0.729; %1         % Inertia Weight
wdamp=1;%PSOdata.wdamp;          % Inertia Weight Damping Ratio
c1=0;%1;%PSOdata.c1;         % Personal Learning Coefficient
c2=0;%3;%PSOdata.c2;             % Global Learning Coefficient

% Velocity Limits
VelMax= PSOdata.VelMax;
VelMin=-VelMax;

nParticleMutation = PSOdata.nParticleMutation;%2 + 3;%3 8;    % Number of Mutations Performed on Each Particle
nGlobalBestMutation = PSOdata.nGlobalBestMutation;%5 + 5;%5 10;  % Number of Mutations Performed on Global Best

%% Initialization
GlobalBest = BestSol;
stuckCounter = 0;

for i=1:nPop
    particle(i).Position = GlobalBest.Position;
end

%% PSO Main Loop
OldMaxIt = it;
NewMax = it + 7500;

for it=OldMaxIt:NewMax
    
    for i=1:nPop        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Velocity Mirror Effect
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        % Evaluation
        [particle(i).Cost, particle(i).Sol] = CostFunction(particle(i).Position);
        
        % Perform Mutation
        for j=1:nParticleMutation
            NewParticle = particle(i);
            NewParticle.Position = PSO_Mutate(particle(i).Position);
            [NewParticle.Cost, NewParticle.Sol] = CostFunction(NewParticle.Position);
            if NewParticle.Cost <= particle(i).Cost
                particle(i) = NewParticle;
            end
        end
        
        % Update Personal Best
        if particle(i).Cost<=particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            particle(i).Best.Sol=particle(i).Sol;
            
            % Update Global Best
            if particle(i).Best.Cost<=GlobalBest.Cost
                GlobalBest=particle(i).Best;
            end
            
        end
        
    end
    
    % Perform Mutation on Global Best
    for i=1:nGlobalBestMutation
        NewParticle = GlobalBest;
        NewParticle.Position = PSO_Mutate(GlobalBest.Position);
        for lm = 1:10
            NewParticle.Position = PSO_Mutate(NewParticle.Position);
            [NewParticle.Cost, NewParticle.Sol] = CostFunction(NewParticle.Position);
            if NewParticle.Cost <= GlobalBest.Cost
                GlobalBest = NewParticle;
            end
        end
        [NewParticle.Cost, NewParticle.Sol] = CostFunction(NewParticle.Position);
        if NewParticle.Cost <= GlobalBest.Cost
            GlobalBest = NewParticle;
        end
    end
    w=w*wdamp;
    
    BestCost(it)=GlobalBest.Cost;
    Costs = [particle.Cost];
    AvgCost(it) = mean(Costs);

    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it)) ...
         ' , Number of Bins = ' num2str(GlobalBest.Sol.nBin)]);      
    
    for i3 = 1:GlobalBest.Sol.nBin
        GlobalBest.Sol.BWeights{i3,1} = model.w(GlobalBest.Sol.BPos{i3}); 
    end

    if (it>=50)  
        if (BestCost(it,1) == BestCost(it-1,1)) && BestCost(it,1) == floor(BestCost(it,1))
            %Stuck in feasible solution
            stuckCounter=stuckCounter+1;
            if (stuckCounter == 100)
                BestIt = it - stuckCounter;
                break;        
            end
            
        elseif (BestCost(it,1) == BestCost(it-1,1))
            %Stuck in unfeasible solution
            stuckCounter=stuckCounter+1;
            if (stuckCounter == 2500)
                BestIt = it - stuckCounter;
                break;        
            end
        else 
            stuckCounter=1;   
        end
        
    end
    
end

BestSol = GlobalBest;

%end
