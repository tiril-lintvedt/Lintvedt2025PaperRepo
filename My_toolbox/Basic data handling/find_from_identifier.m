function index = find_from_identifier(saisir,startpos,str,Option)
%find_from_identifier 		- Finds index of samples with given identifier 
%function [index] = find_from_identifier(saisir,startpos,str,Option)
%gives the indices of saisir rows for which the identifiers contain the string 
% str, in starting position startpos 
%
% Option = -2  % selects from samples based on regular expressions regexp. Replace startpos with a regexp.
% Option = -1  % selects from samples based on the underscore delimiter '_'
% Option = 0  % selects from samples
% Option = 1  % selects from variables


 if (nargin==3) 
        option=0;
 else
        option=Option;
 end    

if (option==0)
    [n,p]=size(saisir.i);
    aux=saisir.i(:,startpos:p);
    index=strmatch(str,aux);
elseif (option==1)
    [n,p]=size(saisir.v);
    aux=saisir.v(:,startpos:p);
    index=strmatch(str,aux);
    
elseif (option == -1)
    name_parts = split(cellstr(saisir.i),'_');
    aux = name_parts(:,startpos); 
    index = strmatch(str,join(aux,'_',2)); % join columns to generalize for multiple name parts input  


elseif (option == -2)
    matches = regexp(cellstr(saisir.i), startpos, 'start');
    index = find(~cellfun(@isempty,matches));     
    
    
end