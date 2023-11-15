function [Solution] = Plotter(BestCost, BestSol)

figure (1)
plot(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');

            
for i = 1:BestSol.Sol.nBin
    BWeights = cell2mat(BestSol.Sol.BWeights(i));
    Solution(1:length(BWeights),i) = BWeights'; %#ok
end

figure (2)
bar(Solution','stacked')
title('Weights assignment')
xlabel('Bins');


end