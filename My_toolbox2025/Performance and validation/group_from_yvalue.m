function [Xsegments,Ysegments,nsampseg] = group_from_value(X,Y,target, nsegments)

% ------------ Divide X block into groups based on Y values --------------
% Divide X into n blocks based on the corresponding values in Y, e.g.
% "Low", "Medium" and "High" reference values. The Value ranges will be of
% equal sizes based on the min-to-max value range. NB! check the resulting  
% number of samples in each group.
% -------------------------------------------------------------------------
%
%   INPUT:
%               X         -  Saisir data structure
%               Y         -  Saisir data structure of corresponding values 
%               target    -  Column number in block Y to group X from 
%                            (In case of multiple Y columns)
%               nsegments  - Number of segments to group X in based on Y
%                            values (i.e number of value ranges)   
%   OUTPUT:
%               Xsegments -  Cell array of saisir data structures
%               nsampseg  -  Vector array remporting how many samples are
%                            in each of the ouput data segments.
% 
% Note:
% *  Assumes row-to-row correspondence between X and Y.
% -------------------------------------------------------------------------
% EXAMPLE CALL
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
    
    a = min(Y.d(:,target));
    b = max(Y.d(:,target)); 
    
    % Defining equal-spaced value ranges (Should have a std based option here as well)
    endPoints = linspace(a,b,nsegments + 1);   %4 endpoints for 3 segments
    startp = endPoints(1:end-1);                %3 starting points
    stopp = endPoints(2:end);                   %3 stopping points
    midPoints = stopp - ((stopp(1)-startp(1))/2); %3 middle points
    
    % Add infinitesimal value to starting point 1. Otherwise at least one 
    % sample will always be excluded in the for loop below
    dy = range([startp(1) stopp(1)])/1000;
    startp(1) = startp(1) - dy;

    nsamseg = [];
    Xsegments = {};
    for i = 1:nsegments
        segind = find((Y.d(:,target) > startp(i)) .* (Y.d(:,target) <= stopp(i)));         
        Xsegments{i} = selectrow(X,segind);
        Ysegments{i} = selectrow(Y,segind);
        nsampseg(i) = length(segind);
    end
    
    
    
end

