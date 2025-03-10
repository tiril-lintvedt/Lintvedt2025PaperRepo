function TM=CreateTransferMatrix(Source,Target)
% TransferredX=SourceX*TM;
% For interpolation (i.e. calibration transfer). 
% Alternative : interp1

Source=Source(:)';
Target=Target(:)';
TM=zeros(numel(Source),numel(Target));

for N=1:numel(Target)
    D=Source-Target(N);
    [MD,ID]=min(abs(D));
    if MD<0.1
        TM(ID,N)=1;
    else
        Before=find(Source<Target(N),1,'last');
        After=find(Source>Target(N),1,'first');
        if ~isempty(Before) && ~isempty(After)
            BD=abs(Source(Before)-Target(N));
            AD=abs(Source(After)-Target(N));
            TM(Before,N)=AD/(BD+AD);
            TM(After,N)=BD/(BD+AD);
        end
    end
    
end