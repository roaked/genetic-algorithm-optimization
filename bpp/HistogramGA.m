%Run multiple times GA

clear all

% Simulation parameters
model = CreateModel(2);
NSimul = 10;

% Initialization
BCosts = zeros(NSimul,1)';
Its = zeros(NSimul,1)';
k = 1;
NN = zeros(2,1)';
time = zeros(NSimul,1)';

% Main loop
for i = 1:NSimul
    tic
    [~,~, BestSol,~, it, GAData] = GA(model);
        time(i) = toc;
    disp(['Simulation ' num2str(i) ': Best Cost = ' num2str(BestSol.Cost) ...
        ': Time = ' num2str(time(i))]);

    BCosts(i) = BestSol.Cost;
    
    if floor(BestSol.Cost) == BestSol.Cost
        NN(k) = BestSol.Cost;
        k = k + 1;
    end
    Its(i) = it;
end

%% Results

if range(NN) == 0
    binc = NN(1);
    counts = length(NN);
else
    binc = min(NN):1:max(NN);
    counts = hist(NN,binc);
end

% Plot of Bin Incidences
figure
X = categorical(binc);
Y = i-k+1;
bar(X,counts)
hold on

if Y > 0
    X1 = categorical({'Non feasible'});
    bar(X1,Y)
end
ylim([0 max(counts)+1])
xlabel('Number of Bins')
ylabel('Number of incidences')
title('Simulations with GA')


%Plot of Best Costs
figure
scatter(1:length(BCosts),BCosts)
xlim([0.5 length(BCosts)+0.5])
%xticks(0:1:length(BCosts))
ylim([min(BCosts)-1 max(BCosts)+1])
hold on

xc = [0 length(BCosts)+1];
yc = [1 1] * model.m;
plot(xc,yc,'--','LineWidth',1);
xlabel('Simulations')
ylabel('Best Costs')
title('Simulations with GA')

%Plot of Number of Iterations
figure
scatter(1:length(Its),Its)
xlim([0.5 length(Its)+0.5])
%xticks(0:1:length(Its))
ylim([0 max(Its)+50])
xlabel('Simulations')
ylabel('Number of Iterations')
title('Simulations with GA')

%Plot of Simulation Time
figure
scatter(1:length(time),time)
xlim([0.5 length(time)+0.5])
%xticks(0:1:length(time))
ylim([0 max(time)+20])
xlabel('Simulations')
ylabel('Time [seconds]')
title('Simulations with GA')