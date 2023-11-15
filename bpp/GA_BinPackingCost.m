
function [J, sol] = GA_BinPackingCost(pos, model)

    n = model.n;
    w = model.w;
    c = model.c;
    
    Sep = find(pos>n);
    
    From = [0 Sep] + 1;
    To = [Sep length(pos)+1] - 1;
    
    B = {};
    for i=1:length(From)
        Bi = pos(From(i):To(i));
        if numel(Bi)>0
            B = [B; Bi]; %#ok
        end
    end
    
    nBin = numel(B);
    Viol = zeros(nBin,1);
    for i=1:nBin
        Vi = sum(w(B{i}));
        Viol(i) = max(Vi/c-1, 0);
    end
    
    MeanViol = mean(Viol);
    
    if n<100
        alpha = 0.7*n;
    else
        alpha = 1.25*n;
    end
    %sum(w)/c = lower bound of bins needed
    J = nBin + alpha*MeanViol + max(sum(model.w)/model.c - nBin,0)*10;
    
    sol.nBin = nBin;
    sol.BPos = B;
    sol.Viol = Viol;
    sol.MeanViol = MeanViol;
end
