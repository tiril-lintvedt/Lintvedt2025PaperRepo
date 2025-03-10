function [C, iAB_C, iAB_Cprime, idx1, idx2] = repintersect(A2, B2)

% Join data for faster sorting:
A = A2 * [4294967296; 1];
B = B2 * [4294967296; 1];
[uA, ~, iAx] = unique(A, "stable");
[a, idx] = sort(iAx);
aV = 1:numel(A);
aV = aV(idx).';
aI = RunLen(a);
[uB,~,iBx] = unique(B, "stable");
[b, idx] = sort(iBx);
bV = 1:numel(B);
bV = bV(idx).';
bL = cumsum([1, RunLen(b)]);
[C2, iuA, iuB] = intersect(uA, uB, "stable");

% Split data for the output:
C = [floor(C2 / 4294967296), rem(C2, 4294967296)];
iAB_C = cell(numel(iuA), 1);
a0 = 1;

for ii = 1:numel(iuA)
   a1 = a0 + aI(ii);       % Easier indexing for A
   aa = aV(a0:a1 - 1);
   a0 = a1;
   na = numel(aa);
   
   b0 = bL(iuB(ii));       % Need accumulated RunLength for B
   b1 = bL(iuB(ii) + 1) - 1;
   bb = bV(b0:b1);
   nb = numel(bb);
   
   qa = repmat(aa.', nb, 1);   % Replace MESHGRID
   qb = repmat(bb,   1, na);
   
   iAB_C{ii} = [qa(:), qb(:)];
end

iAB_Cprime = cell2mat(iAB_C);
idx2 = cumsum(cellfun('size', iAB_C, 1));  % Faster than cellfun(@(x) size(x,1), C)
idx1 = [1; idx2(1:end-1) + 1];
end

% Helper function: 
function n = RunLen(x)
d = [true; diff(x(:)) ~= 0];  % TRUE if values change
n = diff(find([d.', true]));  % Indices of changes
end