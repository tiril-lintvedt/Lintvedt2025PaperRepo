function spec = replace_by_line(spec,pos)
% Replaces points in a spectrum by a line going from spec(pos(1)) to
% spec(pos(end)), according to line definition y = ax + b
% -------------------------------------------------------------------------
% Inputs:
% spec  - spectrum vector
% pos   - index position in spectrum which shall be replaced 
% -------------------------------------------------------------------------
% Output:
% spec  - spectrum vector with positions in "pos" replaced by line
% -------------------------------------------------------------------------

x = [pos(1) pos(end)];
y = [spec(pos(1)) spec(pos(end))]; % Spectrum values at given positions
coeff = [[1; 1]  x(:)]\y(:);   % Calculate Parameter Vector
a = coeff(2);
b = coeff(1);


f = @(x) a*x + b ;

new_vals = f(pos);
spec(pos) = new_vals;


end