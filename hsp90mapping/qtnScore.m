function [qtnVector,candidates] = qtnScore(ph2Mat)

[rows cols]=size(ph2Mat);

qtnVector=0;

%first calculate qtnScore (min of ph2 for each candidate position)
for i=1:rows
   
    toAnalyze=ph2Mat(i,:);
    toAnalyze(toAnalyze==-1)=NaN;
    logP=real(-log10(toAnalyze));
    maxP(i)=min(logP);
    
end

qtnVector=maxP;

%now determine true causal variant(s) from QTN scores

%first determine maximum QTN score variant
[maxQTNscore,maxIdx]=max(qtnVector);

%now check for other candidate variants
%candidates=find(qtnVector>(maxQTNscore-1));
if maxQTNscore>-log10(0.2)
    candidates=find(qtnVector>maxQTNscore*0.7);
else
    candidates=[];
end

end