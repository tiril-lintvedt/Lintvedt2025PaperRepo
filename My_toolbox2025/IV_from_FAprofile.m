function  IV = IV_from_FAprofile(FA_profile, DB_selection)
    % ----- Calculation of Iodine value (IV) from fatty acids profile -----
    % Input : 
    %           FA_profile    - saisir struct, fatty acids profile (g/100g extracted fat)
    %           DB_selection  - str expression, include only FAs with #DB as 
    %                             indicated by expr., ex. '> x' or '< x' or '== x' 
    %
    % Output: 
    %           IV - iodine values for each sample in FA profile
    %  --------------------------------------------------------------------
    % Additional notes:
    % IV (grams of I2 /100g fat ) is a measure for degree of unsaturation 
    % (i.e number of double bonds) in a sample.
    % ---------------------------------------------------------------------
    
    if nargin == 1
        DB_selection = '>0'; % Include FAs with any number of double bonds, if selection is not specified (ORGINAL)
    end
    
    
    nsamples = size(FA_profile.d,1);
    nFAs = size(FA_profile.d,2);
    
    % Set nan values to zero. Probably '<0.1' or other text elements:
    FA_profile.d(isnan(FA_profile.d)) = 0;
    
    C_weight = 12.01;    % g/mol
    H_weight = 1.008;    % g/mol
    O_weight = 16.0;     % g/mol
    I2_weight = 126.9*2; % g/mol % Iodine attaches to double bonds in I2 form
    
    IV = zeros(nsamples,1);
    DBcons = zeros(nsamples,nFAs) ; % Concentration of double bonds for each FA , mol/100g 
    
    for i = 1:nFAs
        FA = FA_profile.v(i,:);
        Output = regexp(FA,'\w*(?<Carbon>\d\d)[-:](?<DoubleBonds>\d)\w*','names'); 
        
        if isempty(Output)   % Skip if not an FA 
            continue
        else
            ncarbon = str2num(Output.Carbon);
            nDB = str2num(Output.DoubleBonds);
            nhydrogen = (ncarbon-2)*2 - 2*nDB + 3 + 3;
            noxygen = 2; % Always two O in a methylated FA (FAMe)
            total_molar_mass = (ncarbon+1)*C_weight + ...   % One extra carbon in methylated version of FA
                               nhydrogen*H_weight + ...
                               noxygen*O_weight;
            
            
            % Skip if number of double bonds in FA is not appropriate 
            if eval([num2str(nDB) , DB_selection ]) 
                DBcons(:,i) = (FA_profile.d(:,i).*nDB)/total_molar_mass ; 
                
            else 
                continue 
            end
                          
        end
    end
    
    IV = I2_weight.* sum(DBcons,2);


end

