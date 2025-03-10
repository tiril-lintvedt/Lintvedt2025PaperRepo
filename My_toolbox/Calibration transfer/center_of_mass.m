function CM=center_of_mass(X)

P=1:size(X,2);
for N=1:size(X,1)
    CM(N)=(X(N,:)-mean(X(N,:)))*P'/sum(X(N,:));
end