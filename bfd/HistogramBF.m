%Run multiple times Best Fit

clear all

model = CreateModel(2);

NN = zeros(100,1)';

for i = 1:100
model.w = model.w(randperm(length(model.w)));
[~, nBins,~,~,~] = BestFit(model);
NN(i) = nBins;
end

%Results
binc = min(NN):1:max(NN);
counts = hist(NN,binc);
figure
bar(binc,counts)
xlabel('Number of Bins')
ylabel('Number of incidences')
