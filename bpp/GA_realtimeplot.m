close all
clear

%% Model selection

model = CreateModel(1);

%% Problem Definition

CostFunction = @(x) GA_BinPackingCost(x, model);  % Objective Function

nVar = 2*model.n-1;     % Number of Decision Variables

%% GA Parameters

if model.n < 100
    MaxIt = 10000;            % Maximum Number of Iterations
else
    MaxIt = 100000;
end

nPop=100;                    % Population Size

pc=0.1; %0.4                 % Crossover Percentage
nc=2*round(pc*nPop/2);       % Number of Offsprings (Parents)

pm=0.9; %0.8                 % Mutation Percentage
nm=round(pm*nPop);           % Number of Mutants

beta=75;                     % Selection Pressure

%Save parameters
GAdata.nPop = nPop;
GAdata.pc = pc;
GAdata.pm = pm;
GAdata.beta = beta;

%% Initialization

%Create Empty Structure
empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Sol=[];

%Create Population Matrix (Array)
pop=repmat(empty_individual,nPop,1);
 
%Initialize Population
for i=1:nPop
    %Initialize Position
    pop(i).Position=randperm(nVar);
    
    %Evaluation
    [pop(i).Cost, pop(i).Sol]=CostFunction(pop(i).Position);   
end

% Sort Population
Costs=[pop.Cost];
[Costs, SortOrder]=sort(Costs);
pop=pop(SortOrder);

% Update Best Solution Ever Found
BestSol=pop(1);

% Update Worst Cost
WorstCost=max(Costs);

% Array to Hold Best Cost Values
BestCost=zeros(50,1);
AvgCost = zeros(50,1);

stuckCounter = 0;

%% GA Main Loop

for it=1:MaxIt
     
    % Calculate Selection Probabilities
    P=exp(-beta*Costs/WorstCost);
    P=P/sum(P);
    
    % Crossover
    popc=repmat(empty_individual,nc/2,2);
    
    for k=1:nc/2    
        % Select Parents
        i1=GA_RouletteWheelSelection(P);
        i2=GA_RouletteWheelSelection(P);
        p1=pop(i1);
        p2=pop(i2);
        
        % Apply Crossover
        [popc(k,1).Position, popc(k,2).Position]=GA_PermutationCrossover(p1.Position,p2.Position);
        
        % Evaluate Offsprings
        [popc(k,1).Cost, popc(k,1).Sol]=CostFunction(popc(k,1).Position);
        [popc(k,2).Cost, popc(k,2).Sol]=CostFunction(popc(k,2).Position);
    end
    popc=popc(:);
    
    % Mutation
    popm=repmat(empty_individual,nm,1);
    for k=1:nm    
        % Select Parent Index
        i=randi([1 nPop]);
        
        % Select Parent
        p=pop(i);
        
        % Apply Mutation
        popm(k).Position=GA_PermutationMutate(p.Position);
        
        % Evaluate Mutant
        [popm(k).Cost, popm(k).Sol]=CostFunction(popm(k).Position);
    end
    
    % Merge Population
    pop=[pop
         popc
         popm]; %#ok

    % Sort Population
    Costs=[pop.Cost];
    [Costs, SortOrder]=sort(Costs);
    pop=pop(SortOrder);
    
    % Truancate Extra Memebrs
    pop=pop(1:nPop);
    Costs=Costs(1:nPop);
    
    % Update Best Solution Ever Found
    BestSol=pop(1);
    
    % Update Worst Cost
    WorstCost=max(WorstCost,max(Costs));
    
    % Update Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    AvgCost(it) = mean(Costs);
        
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))...
        ' , Number of Bins = ' num2str(BestSol.Sol.nBin)]);  
    
    for i3 = 1:BestSol.Sol.nBin
        BestSol.Sol.BWeights{i3,1} = model.w(BestSol.Sol.BPos{i3}); 
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
            if (stuckCounter == 1000)
                break;        
            end
        else 
            stuckCounter=1;   
        end
    end    
end

