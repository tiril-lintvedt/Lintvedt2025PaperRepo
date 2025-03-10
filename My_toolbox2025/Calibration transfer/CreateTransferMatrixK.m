function TM=CreateTransferMatrixK(K,S,S2)

V=[];
for N=1:numel(K)
    V=[V (1:S)'.^(N-1)];
end

I=V*K(:);
TM=zeros(S,S2);
for N=1:S
    if I(N)>=1 && I(N)<=S2
        if ceil(I(N))==I(N)
            TM(N,I(N))=1;
        else
            TM(N,floor(I(N)))=ceil(I(N))-I(N);
            TM(N,ceil(I(N)))=I(N)-floor(I(N));
        end
    end
end