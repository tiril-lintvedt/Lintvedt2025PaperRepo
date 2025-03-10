function [X, Xmeta] = extractSPC_KA(spec)

% ----------- Extracts full SPC information --------------------------------
% Extracts the available Ssectral data and corresponding meta data that can 
% be found in the auditlog of the spec structure (from GSImportspec)
% -------------------------------------------------------------------------
%
%   INPUT: 
%           spec    -   Struct, The full spectral SPC file as obtained by  
%                       the command : spec = GSImportspec('Filepath',0)
%                       (spec is a struct with field "auditlog")
%
%   OUTPUT: 
%           X       -   Spectral data in saisir structure
%           Xmeta   -   Meta data (auditlog) in saisir structure
% 
%   NOTE
%           * The auditlog is extracted based on the form of SPC acquired
%             from the MarqMetrix Raman instrument(old). Might need edit
%             for generalizing to other SPC
% -------------------------------------------------------------------------

% Spectral data to saisir structure ---------------------------------------
X.d = cat(1, spec.data); % convert to saisir structure
X.v = num2str(spec(1).xaxis');
id = split(cellstr(char(spec.name)),'.');
X.i = char(id(:,1));


% Extract metadata --------------------------------------------------------
spec_tab = struct2table(spec); 
spec_meta = spec_tab.auditlog;
%meta_names = regexp(spec_meta,'[A-Z](\w*)(\s*)(\w*)(\s*)(\s*)(\w*)(\s*)((\(*\w*\)*))=', 'match');
%meta_vals = regexp(spec_meta,'[A-Z](\w*)(\s*)(\w*)(\s*)(\s*)(\w*)(\s*)((\(*\w*\)*))=', 'split');

% Untangle metadata 
meta = regexp(spec_meta,'\n', 'split');
meta = reshape([meta{:,1}], size(meta{1},2),size(meta,1) );
meta = [cellfun(@(x) strsplit(x,'='),meta,'UniformOutput',false)];
meta = meta(1:end-1, :);  % % Remove empty cell row. HARDCODED for Kaiser. Generalize betetr later..
nmeta = size(meta,1); % Number of meta  variables
flatmeta = reshape([meta{:}],nmeta*2,[])'; % 296 = 148(nrows in meta)*2 (because there is a set on  meta name and meta value = 1x2 cells)
meta_names = flatmeta(:,1:2:size(flatmeta,2))';
meta_names = meta_names(1:end,1);
meta_vals = flatmeta(:,2:2:size(flatmeta,2));

% Place in regular Matlab saisir struct
Xmeta.d = meta_vals;
Xmeta.v = char(meta_names);
Xmeta.i = X.i;

end