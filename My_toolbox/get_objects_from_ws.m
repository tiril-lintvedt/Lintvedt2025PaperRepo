function out = get_objects_from_ws(classname, workspace)
% Finds all object in workspace of the given class
%
%       INPUT: 
%                   classname   -   str, name of class (e.g. double/char)
%                   workspace   -   str, 'caller' (only current workspace) 
%                                        or 'base' (The workspace where 
%                                        original function call comes from)
% -------------------------------------------------------------------------

     baseVariables = evalin(workspace,'whos');
     out = cell(0);
     for i = 1:length(baseVariables)
         if (strcmpi(baseVariables(i).class , classname)) % compare classnames
             out{length(out) +1}  = baseVariables(i).name;
         end  
     end
end