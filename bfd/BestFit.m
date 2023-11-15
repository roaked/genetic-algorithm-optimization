function [Solution, nBins, time, Ocupation, Space_Left] = BestFit(model)

tic
w = model.w;
c = model.c;

%Best Fit
Space_Left = c;
Solution = 0;
aux = 0;
nBins = 1;

for i = 1:length(w)    
    %Checks where the current item fits best
    for j = 1:length(Space_Left)
        if Space_Left(j) < w(i)
            aux(j) = NaN; %#ok %doesn't fit 
        else
            aux(j) = Space_Left(j) - w(i); %#ok
        end
    end
    [Min,Ind] = min(aux);
    
    if isnan(Min)
        %Creates new bin and assigns the item
        nBins = length(Space_Left) + 1;
        Space_Left(nBins) = c - w(i);
        Solution(1,nBins) = w(i); %#ok
    
    else
        %Places item in the bin (Ind)
        Space_Left(Ind) = Space_Left(Ind) - w(i);
        
        %This section just makes sure one item doesn't overwrite another
        s = size(Solution(:,Ind),1);
        if Solution(s,Ind) == 0
            k = 0;
            assigned = 0;
            while assigned == 0
                k = k + 1;
                if Solution(k,Ind) == 0
                    Solution(k,Ind) = w(i); %#ok
                    assigned = 1;
                end
            end
            
        else
            Solution(s+1,Ind) = w(i); %#ok
        end
    end
end

%Performance evaluation
time = toc;
soma = sum(w);
Ocupation = soma / (c * length(Space_Left)) * 100;

end