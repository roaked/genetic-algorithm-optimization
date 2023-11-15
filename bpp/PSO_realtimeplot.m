close all
clear

%% Model selection

model = CreateModel(1);

%% Problem Definition

CostFunction = @(x) PSO_BinPackingCost(x, model);  % Objective Function

nVar = 2*model.n-1;     % Number of Decision Variables
VarSize = [1 nVar];     % Decision Variables Matrix Size

VarMin = 0;     % Lower Bound of Decision Variables
VarMax = nVar;     % Upper Bound of Decision Variables


%% PSO Parameters

if model.n < 100
    MaxIt = 5000;            % Maximum Number of Iterations
else
    MaxIt = 75000;
end

nPop=10;            % Population Size (Swarm Size)

w=0.9;%0.729; %1      % Inertia Weight
wdamp=0.99;          % Inertia Weight Damping Ratio
c1=2.0;         % Personal Learning Coefficient
c2=2.0;             % Global Learning Coefficient

% Velocity Limits
VelMax=0.0001; %default 0.02
VelMin=-VelMax;

nParticleMutation = 5; %2 + 8;    %3;    % Number of Mutations Performed on Each Particle
nGlobalBestMutation = 10; %5 + 10; %5;  % Number of Mutations Performed on Global Best

%Save parameters
PSOdata.nPop = nPop;
PSOdata.w = w;
PSOdata.wdamp = wdamp;
PSOdata.c1 = c1;
PSOdata.c2 = c2;
PSOdata.VelMax = VelMax;
PSOdata.nParticleMutation = nParticleMutation;
PSOdata.nGlobalBestMutation = nGlobalBestMutation;

%% Initialization

empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Sol=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
empty_particle.Best.Sol=[];

particle=repmat(empty_particle,nPop,1);

GlobalBest.Cost=inf;

for i=1:nPop
    % Initialize Position
    particle(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    % Initialize Velocity
    particle(i).Velocity=zeros(VarSize);
    
    % Evaluation
    [particle(i).Cost, particle(i).Sol]=CostFunction(particle(i).Position);
    
    % Update Personal Best
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    particle(i).Best.Sol=particle(i).Sol;
    
    % Update Global Best
    if particle(i).Best.Cost<GlobalBest.Cost
        GlobalBest=particle(i).Best;
    end
end

BestCost=zeros(50,1);
AvgCost = zeros(50,1);
stuckCounter = 0;

%% PSO Main Loop

for it=1:MaxIt
    
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
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            particle(i).Best.Sol=particle(i).Sol;
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                GlobalBest=particle(i).Best;
            end
            
        end
        
    end
    
    % Perform Mutation on Global Best
    for i=1:nGlobalBestMutation
        NewParticle = GlobalBest;
        NewParticle.Position = PSO_Mutate(GlobalBest.Position);
        [NewParticle.Cost, NewParticle.Sol] = CostFunction(NewParticle.Position);
        if NewParticle.Cost <= GlobalBest.Cost
            GlobalBest = NewParticle;
        end
    end
    w=w*wdamp;
    
    BestSol = GlobalBest; 
    BestCost(it) = GlobalBest.Cost;
    Costs = [particle.Cost];
    AvgCost(it) = mean(Costs);

    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it)) ...
         ' , Number of Bins = ' num2str(GlobalBest.Sol.nBin)]);      
    
    for i3 = 1:GlobalBest.Sol.nBin
        BestSol.Sol.BWeights{i3,1} = model.w(GlobalBest.Sol.BPos{i3}); 
    end

    
    Solution = 0;
    for i = 1:BestSol.Sol.nBin
        BWeights = cell2mat(BestSol.Sol.BWeights(i));
        Solution(1:length(BWeights),i) = BWeights'; %#ok
    end
    
    if it == 1 || BestCost(it,1) ~= BestCost(it-1,1)
        figure (1)
        bar(Solution','stacked')
        title('Weights assignment')
        xlabel('Bins');
    end
    
    
    if (it>=50)  
        if (BestCost(it,1) == BestCost(it-1,1)) && BestCost(it,1) == floor(BestCost(it,1))
            %Stuck in feasible solution
            stuckCounter=stuckCounter+1;
            if (stuckCounter == 100)
                break;        
            end
            
        elseif (BestCost(it,1) == BestCost(it-1,1))
            %Stuck in unfeasible solution
            stuckCounter=stuckCounter+1;
            if (stuckCounter == 7500)
                break;        
            end
        else 
            stuckCounter=1;   
        end
        
    end
    
end
