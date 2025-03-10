function  X_scaled = saisir_scale(X, rowscale)

% Scale each row of X according to the given scaling parameter rowscale ------------

[nrow ncol]=size(X.d);
X_scaled.d=zeros(nrow,ncol);

for row= 1:nrow
   if(mod(row,10000)==0)
       disp(num2str([row nrow]));
   end
   X_scaled.d(row,:)=X.d(row,:)/rowscale(row);

end
X_scaled.v=X.v;
X_scaled.i=X.i;


end