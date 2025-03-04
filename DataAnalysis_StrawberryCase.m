% =========================================================================
% PAPER IV - NIR vs Raman robustness comparison
% =========================================================================
%
% ------------------------------------------------------------------------
% 
% Toolboxes required: emsc, GSTools, My_toolbox
%
%
% Author: Tiril Lintvedt
% Nofima, Raw materials and process optimisation
% email address: tiril.lintvedt@nofima.no
% Website: 
% Sept 2023; 
% MATLAB version: R2023a
% OS: Windows 10.0

%------------- BEGIN CODE --------------

%% SETTINGS 
% -------------------------------------------------------------------------

% Main plot colors

plot_colors =[0 0.5 0.9   ; 
              0.9 0.5 0   ; 
              0.2 0.5 0.2 ; 
              0.9 0.3 0.1 ;
              0.7 0.7 0.2 ; 
              0   0   0   ;
              0.5 0.5 0.5 ;
              0.2 0.7 0.7 ;
              0.1 0.3 0.9 ;
              0.3 0.9 0.3 ;
              0.9 0.9 0.9 ;
              0.1 0.1 0.5 ;];

set(0,'defaultfigurecolor',[1 1 1])

%% DATA IMPORT
% -------------------------------------------------------------------------

% Main set data -----------------------------------------------------------
[Xkaiser2021, Xka21month] = load_unscrmat_2_saisir('C:\Users\Tiril Aurora\Documents\01 Nofima\DATA\20231003 NIR vs Raman robustness\Strawberries\Kaiser 2021.mat');
[Xkaiser2022, Xka22month] = load_unscrmat_2_saisir('C:\Users\Tiril Aurora\Documents\01 Nofima\DATA\20231003 NIR vs Raman robustness\Strawberries\Kaiser 2022.mat');
[Xmnir2021, Xmn21month] = load_unscrmat_2_saisir('C:\Users\Tiril Aurora\Documents\01 Nofima\DATA\20231003 NIR vs Raman robustness\Strawberries\MicroNIR 2021.mat');
[Xmnir2022, Xmn22month] = load_unscrmat_2_saisir('C:\Users\Tiril Aurora\Documents\01 Nofima\DATA\20231003 NIR vs Raman robustness\Strawberries\MicroNIR 2022.mat');

% Generalize spectrum ids to match each other and reference names
Xkaiser2021.i = [Xkaiser2021.i, repelem('_',size(Xkaiser2021.d,1),1), num2str(Xka21month.d)];
Xkaiser2021.i = char(replace(cellstr(Xkaiser2021.i),'S','')); 
Xkaiser2021.i = char(replace(cellstr(Xkaiser2021.i),' ','')); 
% Make sure all id numbers have 3 digits 
cellnames = cellstr(Xkaiser2021.i);
newcellnames = cellfun(@(x) regexprep(x,'^(\d\d_)',['0',x(1:3)]),cellnames, 'UniformOutput',false);
Xkaiser2021.i = char(newcellnames);


Xkaiser2022.i = [Xkaiser2022.i, repelem('_',size(Xkaiser2022.d,1),1), num2str(Xka22month.d)];
Xkaiser2022.i = char(replace(cellstr(Xkaiser2022.i),' ','')); % DOUBLE CHECK THIS IS OK
Xkaiser2022.i = char(replace(cellstr(Xkaiser2022.i),'S','')); % DOUBLE CHECK THIS IS OK
% Make sure all id numbers have 3 digits 
cellnames = cellstr(Xkaiser2022.i);
newcellnames = cellfun(@(x) regexprep(x,'^(\d_)',['00',x(1:2)]),cellnames, 'UniformOutput',false);
Xkaiser2022.i = char(newcellnames);
cellnames = cellstr(Xkaiser2022.i);
newcellnames = cellfun(@(x) regexprep(x,'^(\d\d_)',['0',x(1:3)]),cellnames, 'UniformOutput',false);
Xkaiser2022.i = char(newcellnames);

Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'MB','')); % DOUBLE CHECK THIS IS OK
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),' ','')); % DOUBLE CHECK THIS IS OK
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'-','_')); % DOUBLE CHECK THIS IS OK
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'.sam','')); % DOUBLE CHECK THIS IS OK
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'jord_sept21_','')); % DOUBLE CHECK THIS IS OK
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'.','')); % DOUBLE CHECK THIS IS OK
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'_A','A')); % DOUBLE CHECK THIS IS OK
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'A','_A')); % DOUBLE CHECK THIS IS OK
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'_B','B')); % DOUBLE CHECK THIS IS OK
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'B','_B')); % DOUBLE CHECK THIS IS OK
Xmnir2021.v = char(replace(cellstr(Xmnir2021.v),',','.'));
cellnames = cellstr(Xmnir2021.i);
newcellnames = cellfun(@(x) regexprep(x,'(_\d)$',''),cellnames, 'UniformOutput',false);
Xmnir2021.i = char(newcellnames);
Xmnir2021.i = [Xmnir2021.i, repelem('_',size(Xmnir2021.d,1),1), num2str(Xmn21month.d)]; % add month data
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'b','_B'));
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),'a','A'));
Xmnir2021.i = char(replace(cellstr(Xmnir2021.i),' _','_')); % DOUBLE CHECK THIS IS OK
Xmnir2021.i(1171,:) = '096_A_1_202109';
Xmnir2021.i(697,:) =  '017_A_1_202109';


Xmnir2022.i = [Xmnir2022.i, repelem('_',size(Xmnir2022.d,1),1), num2str(Xmn22month.d)];
Xmnir2022.i = char(replace(cellstr(Xmnir2022.i),'jord','')); % DOUBLE CHECK THIS IS OK
Xmnir2022.i = char(replace(cellstr(Xmnir2022.i),' ','')); % DOUBLE CHECK THIS IS OK
Xmnir2022.i = char(replace(cellstr(Xmnir2022.i),'-','_')); % DOUBLE CHECK THIS IS OK
Xmnir2022.i = char(replace(cellstr(Xmnir2022.i),'murano_','')); % DOUBLE CHECK THIS IS OK
Xmnir2022.i = char(replace(cellstr(Xmnir2022.i),'.sam','')); % DOUBLE CHECK THIS IS OK
Xmnir2022.i = char(replace(cellstr(Xmnir2022.i),'A','_A')); % DOUBLE CHECK THIS IS OK
Xmnir2022.i = char(replace(cellstr(Xmnir2022.i),'B','_B')); % DOUBLE CHECK THIS IS OK
Xmnir2022.i = char(replace(cellstr(Xmnir2022.i),'a','_A')); % DOUBLE CHECK THIS IS OK
Xmnir2022.v = char(replace(cellstr(Xmnir2022.v),',','.'));



% Reference data -----------------------------------------------------------

Y = readcell('C:\Users\Tiril Aurora\Documents\01 Nofima\DATA\20231003 NIR vs Raman robustness\Strawberries\Ref.xlsx', 'Sheet', 'Sheet1');
sample_names = Y(2:401,[1,13]); % 14 could be added
nsamp = size(sample_names,1);

References.d = [Y{2:401,10}]'; % Brix only
References.v = char(Y(1,10));
References.i = [char(sample_names(:,1)),repelem('_',nsamp,1), ...
                num2str([sample_names{:,2}]')];  %  ... ,repelem('_',nsamp,1),char(sample_names(:,3))

% Generalize reference id to match spectrum names
References.i = char(replace(cellstr(References.i),' ','')); 
References.i = char(replace(cellstr(References.i),'(2)',''));
References.i = char(replace(cellstr(References.i),'S',''));
cellnames = cellstr(References.i);
newcellnames = cellfun(@(x) regexprep(x,'^(\d\d_)',['0',x(1:3)]),cellnames, 'UniformOutput',false);
References.i = char(newcellnames);
cellnames = cellstr(References.i);
newcellnames = cellfun(@(x) regexprep(x,'^(\d_)',['00',x(1:2)]),cellnames, 'UniformOutput',false);
References.i = char(newcellnames);
cellnames = cellstr(References.i);
newcellnames = cellfun(@(x) regexprep(x,'^(\d\d\d_\d\d\d\d\d\d)$',[x(1:4),'000_',x(5:10)]),cellnames, 'UniformOutput',false);
References.i = char(newcellnames);


% Remove spectra without refernce values from spectral blocks -------------
remove_names = {'041_A_1_202109';
                '041_A_2_202109';
                '041_B_1_202109';
                '041_B_2_202109';
                '041_A_3_202109';
                '041_B_3_202109';
                '058_A_1_202109';
                '058_A_2_202109';
                '058_B_1_202109';
                '058_B_2_202109';
                '058_A_3_202109';
                '058_B_3_202109';
                '065_A_1_202109';
                '065_A_2_202109';
                '065_B_1_202109';
                '065_B_2_202109';
                '065_A_3_202109';
                '065_B_3_202109';
                '060_A_1_202106';
                '060_A_2_202106';
                '060_A_3_202106';
                '060_B_1_202106';
                '060_B_2_202106';
                '060_B_3_202106';
                '066_A_1_202106';
                '066_A_2_202106';
                '066_A_3_202106';
                '066_B_1_202106';
                '066_B_2_202106';
                '066_B_3_202106';
                '043_B_2_202106'; % NIR spectra with nan values (replicate still available)
                '068_A_2_202106' }; % NIR spectra with nan values (replicate still available)

for i = 1:length(remove_names)
    Xkaiser2021 = delete_from_identifier(Xkaiser2021,1,remove_names{i,:});
    Xmnir2021 = delete_from_identifier(Xmnir2021,1,remove_names{i,:});
end

% Rename replicates for consistency across sets
Xkaiser2022.i = char(replace(cellstr(Xkaiser2022.i),'_3_','_1_'));
Xkaiser2022.i = char(replace(cellstr(Xkaiser2022.i),'_4_','_2_'));


%% REMOVE SAMPLES THAT ARE NOT COMMON FOR RAMAN AND NIR SETS
% -------------------------------------------------------------------------

% Remove surplus replicates
Xmnir2021 = delete_from_identifier(Xmnir2021,7,'3');
Xmnir2022 = delete_from_identifier(Xmnir2022,7,'3');


X_subsets_kaiser =  {Xkaiser2021, Xkaiser2022};
X_subsets_mnir = {Xmnir2021, Xmnir2022};

X_subsets_kaiser_intersect = {};
X_subsets_mnir_intersect = {};

for i = 1:length(X_subsets_kaiser)
    
    X_ka_i = X_subsets_kaiser{i};
    X_nir_i = X_subsets_mnir{i};
    
   %  % Check if spec name exist in both sets if not remove spectrum    
   % [common_samples,ia,ib] = intersect(cellstr(X_ka_i.i(:,[1:4, 9:14])),cellstr(X_nir_i.i(:,[1:4, 9:14])));
   %  % Here only the first entry of a match in an array is given in index.
    
    keep_in_kaiser = find(ismember(cellstr(X_ka_i.i(:,[1:4, 9:14])),cellstr(X_nir_i.i(:,[1:4, 9:14]))));
    keep_in_mnir = find(ismember(cellstr(X_nir_i.i(:,[1:4, 9:14])), cellstr(X_ka_i.i(:,[1:4, 9:14]))));
    
    X_subsets_kaiser_intersect{i} = selectrow(X_ka_i,keep_in_kaiser); %#ok<*SAGROW>
    X_subsets_mnir_intersect{i} = selectrow(X_nir_i,keep_in_mnir);

    % idx = 1;
    % for j = 1:length(common_samples)
    %     name = common_samples{j};
    %     id = find_from_identifier(X_ka_i,[1,4],name,-1); 
    %     X_subsets_kaiser_intersect{i}.d(end:(end+length(id)),:) = X_ka_i.d(id,:);
    %     X_subsets_kaiser_intersect{i}.i(end:(end+length(id)),:) = X_ka_i.i(id,:);
    %     X_subsets_kaiser_intersect{i}.v = X_ka_i.v;
    %     idx = idx + length(id);
    % 
    %     % ALSO FOR MNIR......
    %     % X_subsets_mnir_intersect{i} =;   
    % end



    
end

Xkaiser2021 =  X_subsets_kaiser_intersect{1};
Xkaiser2022 = X_subsets_kaiser_intersect{2};
Xmnir2021 = X_subsets_mnir_intersect{1};
Xmnir2022 = X_subsets_mnir_intersect{2};


%% EXPLORE RAW SPECTRA
% -------------------------------------------------------------------------

plot_Raman(str2num(Xkaiser2021.v), Xkaiser2021.d);
plot_Raman(str2num(Xkaiser2022.v), Xkaiser2022.d);
plot_NIR(str2num(Xmnir2021.v), Xmnir2021.d) ;
plot_NIR(str2num(Xmnir2022.v), Xmnir2022.d) ;

%% PREPROCESSING - PART I
% -------------------------------------------------------------------------

% Kaiser 2021 -------------------------------------------------------------

% Raman shift range
[~,i1] = min(abs(str2num(Xkaiser2021.v)-400));  % 520
[~,i2] = min(abs(str2num(Xkaiser2021.v)-1800)); % 1800
Xkaiser2021 = selectcol(Xkaiser2021,i1:i2);

% Kaiser 2022 -------------------------------------------------------------

% Raman shift range
[~,i1] = min(abs(str2num(Xkaiser2022.v)-400));  % 520
[~,i2] = min(abs(str2num(Xkaiser2022.v)-1800)); % 1800
Xkaiser2022 = selectcol(Xkaiser2022,i1:i2);

% MikroNIR ----------------------------------------------------------------
% Keep full region

clear i1 i2

%% EXPLORE NEW SPECTRAL REGIONS
% -------------------------------------------------------------------------

plot_Raman(str2num(Xkaiser2021.v), Xkaiser2021.d);
plot_Raman(str2num(Xkaiser2022.v), Xkaiser2022.d);
plot_NIR(str2num(Xmnir2021.v), Xmnir2021.d) ;
plot_NIR(str2num(Xmnir2022.v), Xmnir2022.d) ;


%% SPIKE REMOVAL FOR RAMAN SPECTRA 
%  ------------------------------------------------------------------------
% No spikes in these spectra
th = 20; % threshold for spike detection
[Xkaiser2021, ~] = spikefix_whitaker_multi(Xkaiser2021,2,2,1,th,1); % double check this later

th = 20; % threshold for spike detection
[Xkaiser2022, ~] = spikefix_whitaker_multi(Xkaiser2022,2,2,1,th,1); % double check this later
% Note that the threshold must be much higher for Kaiser compared to MM to not detect the Raman peaks as spikes 

%% EXPLORE UNSMOOTHED SPECTRA
% -------------------------------------------------------------------------

plot_Raman(str2num(Xkaiser2021.v), Xkaiser2021.d);
plot_Raman(str2num(Xkaiser2022.v), Xkaiser2022.d);

%% SMOOTHING RAMAN SPECTRA
% -------------------------------------------------------------------------

% Kaiser 2021 -------------------------------------------------------------
Xkaiser2021 = saisir_derivative(Xkaiser2021,2,15,0);% Savitsky-Golay smoothing, 2nd degree polynomial, wd size 9.
nvar = length(Xkaiser2021.v);
Xkaiser2021 = selectcol(Xkaiser2021,8:(nvar-7));% Remove edge effects (the four edge points on each side)

% Kaiser 2022 -------------------------------------------------------------
Xkaiser2022 = saisir_derivative(Xkaiser2022,2,15,0);% Savitsky-Golay smoothing, 2nd degree polynomial, wd size 9.
nvar = length(Xkaiser2022.v);
Xkaiser2022 = selectcol(Xkaiser2022,8:(nvar-7));% Remove edge effects (the four edge points on each side)

%% EXPLORE SMOOTHED SPECTRA
% -------------------------------------------------------------------------

plot_Raman(str2num(Xkaiser2021.v), Xkaiser2021.d);
plot_Raman(str2num(Xkaiser2022.v), Xkaiser2022.d);

%% INTENSITY CORRECTION OF RAMAN DATA DUE TO SERVICE ADJUSTMENTS 2022
% -------------------------------------------------------------------------
% During Kaiser instrument service a white reference light was measured before and after 
% the service. There was a discrepancy between the two measurements, due to
% an update in their intensity correction. The spectra below are the white
% lights after the built-in intensity correction.

load('C:\Users\Tiril Aurora\Documents\01 Nofima\DATA\20231003 NIR vs Raman robustness\Wavelength dependent Intensity Correction\PHATPROBEWhitelightservice.mat')

whitelight_service.d = PHATPROBEWhitelightservice;
whitelight_service.v = VarLabels;
whitelight_service.i = ObjLabels;

% Find equal intercept point (Raman shift where spectrum intensity should 
% be unchanged befor and after correction)
tol = 1; % How close the intensity must be in counts to "produce equality" below.
iequal = find(abs(whitelight_service.d(1,:)-whitelight_service.d(2,:))  < tol); 
wnequal = whitelight_service.v(iequal,:);

ph = plot_Raman(str2num(whitelight_service.v), whitelight_service.d);
set(ph, 'Linewidth',2)
xline(str2num(wnequal))
legend(cellstr(whitelight_service.i),'FontSize',14)

% Clip out correct Raman shift region (Corresponding to spectral blocks)
[~,i1] = min(abs(str2num(whitelight_service.v)-360));  % 520
[~,i2] = min(abs(str2num(whitelight_service.v)-1810)); % 1800
whitelight_service = selectcol(whitelight_service,i1:i2);

% Smooth these to avoid introduction of noise during intensity correction
whitelight_service_smooth = saisir_derivative(whitelight_service,1,19,0);% Savitsky-Golay smoothing, 2nd degree polynomial, wd size 9.
nvar = length(whitelight_service.v);
whitelight_service_smooth = selectcol(whitelight_service_smooth,15:(nvar-14));% Remove edge effects (the four edge points on each side)
% Larger window means more edge effects, must remove more. Maybe look at
% larger area before smoothing, then clip to exact area

ph = plot_Raman(str2num(whitelight_service.v), whitelight_service.d); %#ok<NASGU>
hold on
plot(str2num(whitelight_service_smooth.v), whitelight_service_smooth.d)
%set(ph, 'Linewidth',2)
xline(str2num(wnequal))
legend(cellstr(whitelight_service.i),'FontSize',14)

% Clip out correct Raman shift region for white light (Corresponding to spectral blocks)
[~,i1] = min(abs(str2num(whitelight_service_smooth.v)-404));  % Same as Kaiser2022
[~,i2] = min(abs(str2num(whitelight_service_smooth.v)-1796)); % Same as Kaiser2022
whitelight_service_smooth = selectcol(whitelight_service_smooth,i1:i2);



% Adjust the reference before multiplication with spectra, in order to ensure 
% that correction is 1:1 at the point where intensity should be unchanged (wnequal)
% before/after service (cross point between whitelight references)

% Normalize the whitelight ref by the wavenumber of equal intensity 
whitelight_service_smoothnorm = saisir_normpeak(whitelight_service_smooth, str2num(wnequal));

% Correction of spectra
X2021 = Xkaiser2021;
X2022 = Xkaiser2022;

X2022_corr=zeros(800,1393);
for i=1:size(X2022.d,1) 
    X2022_corr(i,:)=X2022.d(i,:)./whitelight_service_smoothnorm.d(2,:);
end

CorrScaledTest=(whitelight_service_smoothnorm.d(1,:)-1)*1; %1.9; % Subrating 1 makes the wnequal go to zero, and ensures that when we multiply, this point is not altered in intensity. 
CorrScaledTest=CorrScaledTest+1; % Then the 1 ust be added to make the wnequal go to 1, which ensures that spectral intensity is not altered upon multiplication
%CorrScaledTest = whitelight_service_smooth.d(1,:); (Simply using the intensiyt corr was not enough to correct for all wl dependent intensity error)

figure
plot(str2num(whitelight_service_smoothnorm.v), CorrScaledTest)
yline(1)
xline(str2num(wnequal))

for i=1:size(X2022.d,1) 
    X2022_corr(i,:)=X2022_corr(i,:).*CorrScaledTest;
end

X2022corr.d = X2022_corr;
X2022corr.v = X2022.v;
X2022corr.i = X2022.i;


%  Compare raw
figure
ph = plot(str2num(X2022.v), X2022_corr, 'r');
hold on
ph2 = plot(str2num(X2022.v), X2022.d, 'b');
legend([ph(1) ph2(1)], 'Corrected 2022', '2022')
yline(0)
yline(1)
xline(str2num(wnequal))

% Compare after baseline correction
X2021bsl = saisir_als(X2021,4,0.001);
X2022bsl = saisir_als(X2022,4,0.001);
X2022corrbsl = saisir_als(X2022corr,4,0.001);

figure;
hold on
%ph2 = plot(str2num(X2022.v), X2022.d, 'b');
ph2 = plot(str2num(X2021bsl.v), X2021bsl.d, 'b');
ph = plot(str2num(X2022bsl.v), X2022corrbsl.d, 'r');
yline(0)
yline(1)
xline(str2num(wnequal))
legend([ph(1) ph2(1)], 'Corrected 2022', '2021')

% Compare after SNV correction
X2021bslsnv = saisir_snv(X2021bsl);
X2022bslsnv = saisir_snv(X2022bsl);
X2022corrbslsnv = saisir_snv(X2022corrbsl);

figure;
hold on
%ph2 = plot(str2num(X2022.v), X2022.d, 'b');
ph2 = plot(str2num(X2021bslsnv.v), X2021bslsnv.d, 'b');
ph = plot(str2num(X2022bslsnv.v), X2022corrbslsnv.d, 'r');
yline(0)
yline(1)
xline(str2num(wnequal))
legend([ph(1) ph2(1)], 'Corrected 2022', '2021')

% Store intensity corrected spectra
Xkaiser2022.d = X2022corr.d;


%% OTHER STUFF THAT WAS TESTED DURING INTENSITY CORRECTION..
% % % This part can be used to optimize the slope of the correction etc..
% % a = 1; % adjustment factor
% % CorrScaledTest=(whitelight_service_smooth.d(1,:))*a;
% % CorrScaledTest=CorrScaledTest;
% % 
% % figure; plot(CorrScaledTest); hold on; 
% % CorrScaledTest=(whitelight_service_smooth.d(1,:)-1)*a;
% % CorrScaledTest=CorrScaledTest+1;
% 
% % Final test: Maybe snv the white refs as JP did. Then
% % the correction factor has more the desired effect i think. The part of
% % spectrum where Intensity should be same must be at 1.
% %3.4 justeringsfaktor
% % a = 1.5; % adjustment factor
% % ScaledRef = whitelight_service_smooth.d;
% % 
% % CorrScaledTest=(snv(whitelight_service_smooth.d(1,:)))*a;
% % CorrScaledTest=CorrScaledTest;
% % CorrScaledTest=(whitelight_service_smooth.d(1,:)-1)*a;
% % CorrScaledTest=CorrScaledTest+1;
% % 
% % figure
% % plot(str2num(whitelight_service_smooth.v), snv(whitelight_service_smooth.d(1,:)))
% % yline(1)
% % xline(str2num(wnequal))
% 
% 
% % diff = whitelight_service_smooth.d(1,:)- whitelight_service_smooth.d(2,:);
% % CorrScaledTest = whitelight_service_smooth.d(1,:)+ 7*diff;
% % diff introduces a lot of noise. Maybe snv the white refs as JP did. Then
% % the correction factor has more the desired effect i think. The part of
% % spectrum where Intensity should be same must be at 1.
% 
% % p = polyfit(1:length(diff),diff, 2)
% % diffline = polyval(p,1:length(diff))
% % hold on; plot(diffline)
% 
% % Adjust the slope with correction factor
% % Make sure the region of equall intensity is still equal. Plot before and
% % after with baseline correction and SNV
% %
% %Xkaiser2022 = X2022_corr;
% 
% % % Regression approach instead?
% % [K,Ri] = MLR(whitelight_service_smooth.d(2,:)',whitelight_service_smooth.d(1,:)');
% % Xwhiteref_corr =whitelight_service_smooth.d(2,:)*K(2:end)+K(1);
% % figure;
% % hold on
% % plot(whitelight_service_smooth.d(2,:),'-')
% % plot(Xwhiteref_corr,'--')
% % plot(whitelight_service_smooth.d(1,:), '-');
% % legend('after', 'after transformed to before','before')
% % X2021 = Xkaiser2021;
% % X2022 = Xkaiser2022;
% % X2022_corr.d = X2022.d*K(2:end)+K(1);
% % X2022_corr.v = X2022.v;
% % X2022_corr.i = X2022.i;
% % figure;
% % hold on
% % p1 = plot(str2num(X2022.v),X2022.d,'b-');
% % p2 = plot(str2num(X2022.v), X2022_corr.d,'r-');
% % legend([p1(1) p2(1)],'2022', '2022 corrected')
% % 
% % 
% % % Check some prep
% % X2022_als = saisir_als(X2022,4,0.001);
% % X2022_corr_als = saisir_als(X2022_corr,4,0.001);
% % X2021_als = saisir_als(X2021,4,0.001);
% % 
% % % Doesn't seem to work well enough. NOt surprising wrt the complexity of
% % % the MLR model.. doesn't incorporate wavelrngth specific correction..
% 
% % % Transformation matrix approach instead ?
% T_after2before = whitelight_service_smooth.d(1,i1:i2)'/whitelight_service_smooth.d(2,i1:i2)';
% WhitelightAfter_corr = T_after2before*whitelight_service_smooth.d(2,i1:i2)';
% %whitelight_service_smooth.d(1,i1:i2)'*inv(A)
% figure;
% hold on
% plot(whitelight_service_smooth.d(1,i1:i2),'-')
% plot(whitelight_service_smooth.d(2,i1:i2),'-')
% plot(WhitelightAfter_corr, '--');
% legend('before', 'after','after transformed to before')
% 
% %
% X2021 = Xkaiser2021;
% X2022 = Xkaiser2022;
% X2022_corr_1 = T_after2before*X2022.d(1,:)';
% 
% Ta2b_diag = diag(T_after2before(:,1393));
% X2022_corr_1 = Ta2b_diag*X2022.d(1,:)';
% figure;plot(X2022_corr_1);hold on
% plot(X2022.d(1,:))
% 
% 
% figure
% hold on
% p1 = plot(str2num(X2021.v), X2021.d,'color', plot_colors(1,:));
% p2 = plot(str2num(X2022.v), X2022_corr,'color', plot_colors(2,:));
% p3 = plot(str2num(X2022.v), X2022.d,'color', plot_colors(3,:));
% legend([p1, p2],'2021 original', '2022 corrected')

%% SEPARATE INSTRUMENTS AND SUBSETS
% -------------------------------------------------------------------------

% Different region versions
X_instrument_subsets = {Xkaiser2021, Xkaiser2022, Xmnir2021, Xmnir2022};
X_subsets = {};
n_subsets = {};

for i = 1:size(X_instrument_subsets,2)

    X = X_instrument_subsets{i};
    % Separate by muscle/part
    [X_sunnyside, ~] = select_from_identifier(X,5,'A'); 
    [X_shadowside, ~] = select_from_identifier(X,5,'B'); 

    % Gather all data subsets ------------------------------------------------- 
    nsets = 3;
    X_subsets(end+1:end+nsets) = {X_sunnyside, X_shadowside, X};
    n_subsets(end+1 : end+nsets , 1) = {'Sunny side'; 'Shadow side'; 'All'};

end



keep X_subsets n_subsets plot_colors References

%% CONNECT REFERENCE DATA WITH RAMAN SPECTRA
% -------------------------------------------------------------------------

ref_idpos = [1:4,9:14]; % string position in reference id tags
x_idpos = [1:4,9:14]; %string position in spectral/data blocks 
reference_names = References.i(:,ref_idpos);
nref = size(References.d,2); % Number of reference analyses
Y_subsets = {};

for i_set = 1:length(X_subsets)  % For each subset, create an R block

    Xi = X_subsets{i_set};
    R.d = zeros(size(Xi.i,1),nref);
    R.v = References.v;
    R.i = {};

    for i_sample = 1:size(Xi.i,1)
        name = Xi.i(i_sample,x_idpos);

        for ref = 1:size(reference_names,1)
            ref_name = reference_names(ref,:);       
            if strcmp(name,ref_name)
                R.d(i_sample,:) = References.d(ref,:);
                R.i{i_sample} = Xi.i(i_sample,:);
                break
            end     
        end

    end
    R.i = char(R.i');
    Y_subsets{i_set} = R;

    % Check if all spectra has a reference value assigned in Y block
    nanflag = any(isnan(Y_subsets{i_set}.d(:))); 
    if nanflag 
        warning(['NAN values were detected in Y_subsets{',num2str(i_set),'}.'])
    end

end



%% BASELINE CORRECTION ONLY (ALS) FOR RAMAN SPECTRA
% -------------------------------------------------------------------------

X_prep_als_subsets = {};
lambda = 4 ;% Smoothing parameter % 4
p = 0.001;

for i = 1:6 % Kaiser 2021 + 2022
    Subset = X_subsets{i};
    [Subset_prep,baseline,wgts] = saisir_als(Subset, lambda, p); % correct baseline
    X_prep_als_subsets{i} = Subset_prep;
    
    % Control plot
    plot_Raman(str2num(Subset_prep.v), Subset_prep.d);

end

%% COMPARE 2021 AND 2022 RAMAN SPECTRA (BEFORE AND AFTER SERVICE)
% -------------------------------------------------------------------------
X2021 = X_prep_als_subsets{3};
X2022 = X_prep_als_subsets{6};

ph = plot_Raman(str2num(X2021.v), X2021.d) ;
set(ph,'color', plot_colors(1,:))
hold on
p22 = plot(str2num(X2022.v), X2022.d, 'Color',plot_colors(2,:));
%plot(str2num(X2022.v), X2022.d);

% Sugar related peaks:
% Look at Olgas Apple paper supplementary
xline(1078)
xline(1124)
xline(918)
xline(629)
xline(493)
xline(834)

leg = legend([ph(1),p22(1)], '2021','2022', 'FontSize', 14, 'Location', 'Northwest');
title(leg, 'Year', 'FontSize', 14)

%% BASELINE CORRECTION (ALS) and SNV FOR RAMAN SPECTRA
% -------------------------------------------------------------------------

% Kaiser ------------------------------------------------------------------

X_prep_alssnv_subsets = {};
lambda = 4; % Smoothing parameter % 4
p = 0.001;

for i = 1:6 % Kaiser 2021 + 2022
    Subset = X_subsets{i};
    [Subset_prep,baseline,wgts] = saisir_als(Subset, lambda, p); % correct baseline
    Subset_prep = saisir_snv(Subset_prep);
    X_prep_alssnv_subsets{i} = Subset_prep;
    
    % Control plot
    plot_Raman(str2num(Subset_prep.v), Subset_prep.d);
    color_by(Y_subsets{i}.d)

end


% MikroNIR ------------------------------------------------------------------

lambda = 1; % Smoothing parameter % 4
p = 0.001;

for i = 7:12 % Kaiser 2021 + 2022
    Subset = X_subsets{i};
    [Subset_prep,baseline,wgts] = saisir_als(Subset, lambda, p); % correct baseline
    %Subset_prep = saisir_snv(Subset_prep);
    %Subset_prep = saisir_normpeak(Subset_prep, 963);
    X_prep_alssnv_subsets{i} = Subset_prep;
    
    % Control plot
    plot_NIR(str2num(Subset_prep.v), Subset_prep.d);
    color_by(Y_subsets{i}.d);

end

%% COMPARE 2021 AND 2022 RAMAN SPECTRA (BEFORE AND AFTER SERVICE)
% -------------------------------------------------------------------------
X2021 = X_prep_alssnv_subsets{3};
X2022 = X_prep_alssnv_subsets{6};

ph = plot_Raman(str2num(X2021.v), X2021.d) ;
set(ph,'color', plot_colors(1,:))
hold on
p22 = plot(str2num(X2022.v), X2022.d, 'Color',plot_colors(2,:));

% Sugar related peaks:
% Look at Olgas Apple paper supplementary
xline(1458)
xline(1265)
xline(1078)
xline(1124)
xline(918)
xline(629)
xline(493)
xline(834)


leg = legend([ph(1),p22(1)], '2021','2022', 'FontSize', 14, 'Location', 'Northwest');
title(leg, 'Year', 'FontSize', 14)

%% BASELINE CORRECTION (ALS) AND NORM. BY INTERNAL STANDARD FOR RAMAN SPECTRA
% -------------------------------------------------------------------------

% Kaiser ------------------------------------------------------------------

X_prep_alsnormpeak_subsets = {};
lambda = 4; % Smoothing parameter % 4
p = 0.001;
wn = 1556 ; % 605 (proximately same as doing SNV);

for i = 1:6 % Kaiser 2021 + 2022
    Subset = X_subsets{i};
    [Subset_prep,baseline,wgts] = saisir_als(Subset, lambda, p); % correct baseline
    Subset_prep = saisir_normpeak(Subset_prep,wn);
    X_prep_alsnormpeak_subsets{i} = Subset_prep;
    
    % Control plot
    plot_Raman(str2num(Subset_prep.v), Subset_prep.d);
    color_by(Y_subsets{i}.d)

end


%% COMPARE 2021 AND 2022 RAMAN SPECTRA (BEFORE AND AFTER SERVICE)
% -------------------------------------------------------------------------
X2021 = X_prep_alsnormpeak_subsets{3};
X2022 = X_prep_alsnormpeak_subsets{6};

ph = plot_Raman(str2num(X2021.v), X2021.d) ;
set(ph,'color', plot_colors(1,:))
hold on
p22 = plot(str2num(X2022.v), X2022.d, 'Color',plot_colors(2,:));

% Sugar related peaks:
% Look at Olgas Apple paper supplementary
xline(1458)
xline(1265)
xline(1078)
xline(1124)
xline(918)
xline(629)
xline(493)
xline(834)


leg = legend([ph(1),p22(1)], '2021','2022', 'FontSize', 14, 'Location', 'Northwest');
title(leg, 'Year', 'FontSize', 14)


%% SNV FOR NIR SPECTRA
% -------------------------------------------------------------------------

X_prep_snv_subsets = {};

for i = 7:12
    Subset = X_subsets{i};
    [Subset_prep] = saisir_snv(Subset); % SNV correction
    X_prep_snv_subsets{i} = Subset_prep;
    
    % Control plot
    plot_NIR(str2num(Subset_prep.v), Subset_prep.d);
    color_by(Y_subsets{i}.d(:,1));

end

%% DERIVATIVE FOR NIR (2nd)
% -------------------------------------------------------------------------

X_prep_deriv2_subsets = {};
derorder = 2; % Removes baseline and multiplicative effects
wdsz = 9;

for i = 7:12
    Subset = X_subsets{i};
    [Subset_prep] = saisir_derivative(Subset,2,wdsz,derorder); % SG 
    Subset_prep =  selectcol(Subset_prep, wdsz:(size(Subset_prep.d,2)-wdsz)) ; % Remove lost endpoints 
    X_prep_deriv2_subsets{i} = Subset_prep;
    
    % Control plot
    plot_NIR(str2num(Subset_prep.v), Subset_prep.d);
    color_by(Y_subsets{i}.d(:,1));

end

%% DERIVATIVE FOR NIR (1ST)
% -------------------------------------------------------------------------

X_prep_deriv1_subsets = {};
derorder = 1; % Removes baseline and multiplicative effects
wdsz = 9;

for i = 7:12
    Subset = X_subsets{i};
    [Subset_prep] = saisir_derivative(Subset,2,wdsz,derorder); % SG 
    Subset_prep =  selectcol(Subset_prep, wdsz:(size(Subset_prep.d,2)-wdsz)) ; % Remove lost endpoints 
    X_prep_deriv1_subsets{i} = Subset_prep;
    
    % Control plot
    plot_NIR(str2num(Subset_prep.v), Subset_prep.d);
    color_by(Y_subsets{i}.d(:,1));

end

%% DERIVATIVE FOR NIR (2ND) + SNV
% -------------------------------------------------------------------------

X_prep_deriv2snv_subsets = {};
derorder = 2; % Removes baseline and multiplicative effects
wdsz = 9;

for i = 7:12
    Subset = X_subsets{i};
    [Subset_prep] = saisir_derivative(Subset,2,wdsz,derorder); % SG 
    Subset_prep =  selectcol(Subset_prep, wdsz:(size(Subset_prep.d,2)-wdsz)) ; % Remove lost endpoints 
    Subset_prep = saisir_snv(Subset_prep);
    X_prep_deriv2snv_subsets{i} = Subset_prep;
   
    
    % Control plot
    plot_NIR(str2num(Subset_prep.v), Subset_prep.d);
    color_by(Y_subsets{i}.d(:,1));
 

end

%% DERIVATIVE FOR NIR (1ST) + SNV
% -------------------------------------------------------------------------

X_prep_deriv1snv_subsets = {};
derorder = 1; % Removes baseline and multiplicative effects
wdsz = 9;

for i = 7:12
    Subset = X_subsets{i};
    [Subset_prep] = saisir_derivative(Subset,2,wdsz,derorder); % SG 
    Subset_prep =  selectcol(Subset_prep, wdsz:(size(Subset_prep.d,2)-wdsz)) ; % Remove lost endpoints 
    Subset_prep = saisir_snv(Subset_prep);
    X_prep_deriv1snv_subsets{i} = Subset_prep;
   
    
    % Control plot
    plot_NIR(str2num(Subset_prep.v), Subset_prep.d);
    color_by(Y_subsets{i}.d(:,1));
 

end

%% COMPARE 2021 AND 2022 MIKRONIR SPECTRA
% -------------------------------------------------------------------------

% Before SNV --------------------------------------------------------------

X2021 = X_subsets{9};
X2022 = X_subsets{12};

ph = plot_NIR(str2num(X2021.v), X2021.d) ;
set(ph,'color', plot_colors(1,:))
hold on
p22 = plot(str2num(X2022.v), X2022.d, 'Color',plot_colors(2,:));
leg = legend([ph(1),p22(1)], '2021','2022', 'FontSize', 14, 'Location', 'Northwest');
title(leg, 'Year', 'FontSize', 14)



% After SNV ---------------------------------------------------------------

X2021 = X_prep_snv_subsets{9};
X2022 = X_prep_snv_subsets{12};

ph = plot_NIR(str2num(X2021.v), X2021.d) ;
set(ph,'color', plot_colors(1,:))
hold on
p22 = plot(str2num(X2022.v), X2022.d, 'Color',plot_colors(2,:));
leg = legend([ph(1),p22(1)], '2021','2022', 'FontSize', 14, 'Location', 'Northwest');
title(leg, 'Year', 'FontSize', 14)
title('SNV')

%% SPECTRA PER MONTH - AFTER PREPROCESSING
% -------------------------------------------------------------------------

X_snv = saisir_rowmerge(X_prep_snv_subsets{9},X_prep_snv_subsets{12});

plot_from_identifier(X_snv,9:14)
ylabel('Normalized absorbance', 'FontSize',18)
xlabel('Wavelength (nm)', 'FontSize',18)
grid on
set(gcf,'renderer','Painters')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/NIRspectra_month_prep','png')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/NIRspectra_month_prep','epsc')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/NIRspectra_month_prep','fig')

X_snv_avg = average_from_identifier(X_snv,9:14);
plot_from_identifier(X_snv_avg,1:6)
ylabel('Normalized absorbance', 'FontSize',18)
xlabel('Wavelength (nm)', 'FontSize',18)
grid on
set(gcf,'renderer','Painters')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/NIRspectraAvg_month_prep','png')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/NIRspectraAvg_month_prep','epsc')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/NIRspectraAvg_month_prep','fig')

% RAMAN

X_snv = saisir_rowmerge(X_prep_alssnv_subsets{3},X_prep_alssnv_subsets{6});

plot_from_identifier(X_snv,9:14)
ylabel('Normalized intensity', 'FontSize',18)
xlabel('Raman shift (cm^{-1})', 'FontSize',18)
grid on

set(gcf,'renderer','Painters')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSspectra_month_prep','png')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSspectra_month_prep','epsc')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSspectra_month_prep','fig')

X_snv_avg = average_from_identifier(X_snv,9:14);
plot_from_identifier(X_snv_avg,1:6)
ylabel('Normalized intensity', 'FontSize',18)
xlabel('Raman shift (cm^{-1})', 'FontSize',18)
grid on

set(gcf,'renderer','Painters')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSspectraAvg_month_prep','png')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSspectraAvg_month_prep','epsc')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSspectraAvg_month_prep','fig')



%% SPECTRAL PCA - EXPLORING VARIATIONS IN SUGARS - RAMAN
% -------------------------------------------------------------------------

% Merge 2021 and 2022 data blocks
%X_als = saisir_rowmerge(X_prep_als_subsets{3},X_prep_als_subsets{6}); %#ok<NASGU>
X_alssnv = saisir_rowmerge(X_prep_alssnv_subsets{3},X_prep_alssnv_subsets{6}); %#ok<NASGU>
%X_alsnormpeak = saisir_rowmerge(X_prep_alsnormpeak_subsets{3},X_prep_alsnormpeak_subsets{6}); 
Y = saisir_rowmerge(Y_subsets{3},Y_subsets{6}); %#ok<NASGU>

% % Select a subsets
% X_alsnormpeak_2021 = select_from_identifier(X_alsnormpeak,9, '2021');
% X_alsnormpeak_2022 = select_from_identifier(X_alsnormpeak,9, '2022');
% 
% X_alsnormpeak_202106 = select_from_identifier(X_alssnv,9, '202106');
% X_alssnv_202109 = select_from_identifier(X_alssnv,9, '202109');
% X_alssnv_202206 = select_from_identifier(X_alssnv,9, '202206');
% X_alssnv_202208 = select_from_identifier(X_alssnv,9, '202208');

ncomp = 6;
[T,cum_var,var,mx,V,U] = pca(X_alssnv.d,ncomp);
plot_comp(T,cum_var,var,mx,V,U,ncomp,X_alssnv.v)
plot_scores(T,X_alssnv, 9:14,'colorBy',Y.d) %9:14, %5
set(gcf,'Position',[200 200 800 900])
set(gcf,'renderer','Painters')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/StrawberryRamanPCA','png')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/StrawberryRamanPCA','epsc')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/StrawberryRamanPCA','fig')


% Plot scores by sugar content
figure
hold on
for i = 1:length(T)
    plot(T(i,1)*1, T(i,2),'*')
end
cbh = color_by(Y.d);
title(cbh,'Brix')
xlabel('Score PC1', 'Fontsize', 18)
ylabel('Score PC2', 'Fontsize', 18)


%% SPECTRAL PCA - EXPLORING VARIATIONS IN SUGARS - NIR
% -------------------------------------------------------------------------

% Merge 2021 and 2022 data blocks

X_snv = saisir_rowmerge(X_prep_snv_subsets{9},X_prep_snv_subsets{12});
% X_deriv2 = saisir_rowmerge(X_prep_deriv2_subsets{9},X_prep_deriv2_subsets{12}); 
% X_alssnv = saisir_rowmerge(X_prep_alssnv_subsets{9},X_prep_alssnv_subsets{12});% 
Y = saisir_rowmerge(Y_subsets{9},Y_subsets{12});

% % Select a subsets
% X_snv_2021 = select_from_identifier(X_snv,9, '2021');
% X_snv_2022 = select_from_identifier(X_snv,9, '2022');
% 
% X_snv_202106 = select_from_identifier(X_snv,9, '202106');
% X_snv_202109 = select_from_identifier(X_snv,9, '202109');
% X_snv_202206 = select_from_identifier(X_snv,9, '202206');
% X_snv_202208 = select_from_identifier(X_snv,9, '202208');

ncomp = 6;
[T,cum_var,var,mx,V,U] = pca(X_snv.d,ncomp);
plot_comp(T,cum_var,var,mx,V,U,ncomp,X_snv.v)
plot_scores(T,X_snv, 9:14,'colorBy',Y.d) %9:14, %5
set(gcf,'Position',[200 200 800 900])

set(gcf,'renderer','Painters')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/StrawberryNIRPCA','png')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/StrawberryNIRPCA','epsc')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/StrawberryNIRPCA','fig')



% Plot scores by sugar content
figure
hold on
for i = 1:length(T)
    plot(T(i,1)*1, T(i,2),'*')
end
cbh = color_by(Y.d);
title(cbh,'Brix')
xlabel('Score PC1', 'Fontsize', 18)
ylabel('Score PC2', 'Fontsize', 18)

% Plot scores by sugar content
figure
hold on
for i = 1:length(T)
    plot(T(i,3)*1, T(i,6),'*')
end
cbh = color_by(Y.d);
title(cbh,'Brix')
xlabel('Score PC3', 'Fontsize', 18)
ylabel('Score PC6', 'Fontsize', 18)



%% LEAVE-ONE-OUT CROSS VALIDATION, RAMAN
% -------------------------------------------------------------------------
% To show that Raman and NIR have similar results, so that it is a good
% case to compare the dependence on number of samples on.

% Merge 2021 and 2022 data blocks
%X_als = saisir_rowmerge(X_prep_als_subsets{3},X_prep_als_subsets{6}); 
X_alssnv = saisir_rowmerge(X_prep_alssnv_subsets{3},X_prep_alssnv_subsets{6});
%X_alsnormpeak = saisir_rowmerge(X_prep_alsnormpeak_subsets{3},X_prep_alsnormpeak_subsets{6}); 

Y = saisir_rowmerge(Y_subsets{3},Y_subsets{6});

idpos = [1:3,9:14];
aopt = 20;
target = 1;
[Predictions,Regcoeff, Performance] = PLSR_CV_from_identifier(X_alssnv, Y, target, idpos, aopt,'PlotLimits',[1 20]);

figs = findobj('type','Figure');
FIG1 = figs(2);
FIG2 = figs(1);

saveas(FIG1,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSstrawb_predvstarget_LooCV_rmseVSncomp','png')
saveas(FIG1,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSstrawb_predvstarget_LooCV_rmseVSncomp','epsc')
saveas(FIG1,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSstrawb_predvstarget_LooCV_rmseVSncomp','fig')
saveas(FIG2,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSstrawb_predvstarget_LooCV_predvstarget','png')
saveas(FIG2,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSstrawb_predvstarget_LooCV_predvstarget','epsc')
saveas(FIG2,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSstrawb_predvstarget_LooCV_predvstarget','fig')



plot_regcoef(X_alssnv.v,squeeze(Regcoeff.d(:,:,4)),0.5, 'findPeaks','no')

%% VALIDATION 2021 -> 2022 RAMAN // OR BY MONTH
% -------------------------------------------------------------------------

% Merge 2021 and 2022 data blocks
%X_als = saisir_rowmerge(X_prep_als_subsets{3},X_prep_als_subsets{6}); 
X_alssnv = saisir_rowmerge(X_prep_alssnv_subsets{3},X_prep_alssnv_subsets{6});
%X_alsnormpeak = saisir_rowmerge(X_prep_alsnormpeak_subsets{3},X_prep_alsnormpeak_subsets{6}); 
Y = saisir_rowmerge(Y_subsets{3},Y_subsets{6});

% % Look only a side A/B, one at a time
% side = 'A'; % 'B'
% X_alssnv_A = select_from_identifier(X_alssnv,5,side); 
% Y_A = select_from_identifier(Y,5,side); %#ok<NASGU>

% % Restric region
%  X_alssnv = selectwn(X_alssnv, 440:1500);

% Validate from identifier (year)
aopt = 15;
[FIG, Regcoeff, Performance] = PLSR_validate_from_identifier2(X_alssnv,Y,1,...
                                9:14,aopt,[1:3,9:14], 'aoptAlg','WM', ...
                                'aoptAlgMetric','RMSEP', 'plotLimits', [1 20],...
                                'aoptMode','calCV'); %#ok<ASGLU> %      %9:12(14) for year only  % LV 5-6 Best choice Raman with test set. CHck later if this is different number when decided from crossval on cal set.
lgd = findobj('type','legend');
delete(lgd)
set(gcf,'renderer','Painters')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSstrawb_predvstarget_month','png')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSstrawb_predvstarget_month','epsc')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSstrawb_predvstarget_month','fig')


% Plot performance summary
[FIG3, ah] = plot_allcombos_performance(Performance);
set(ah{1}, 'YLim',[0, 1.7])
set(ah{2}, 'YLim',[-2 2])
set(ah{3}, 'YLim',[0, 1.5])

set(gcf,'renderer','Painters')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSstrawb_errors_month','png')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSstrawb_errors_month','epsc')
saveas(gcf,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSstrawb_errors_month','fig')



plot_regcoef(Regcoeff.v, Regcoeff.d, 0.5,'findPeaks','no','plotColors', plot_colors)
xlabel('Raman shift (cm^{-1})', 'FontSize',18)

% Sugar related peaks:
xline(1458, 'Color', [0.5 0.85 0.5],'Linewidth', 3)
xline(1265, 'Color', [0.5 0.85 0.5],'Linewidth', 3)
xline(1078, 'Color', [0.5 0.85 0.5],'Linewidth', 3)
xline(1124, 'Color', [0.5 0.85 0.5],'Linewidth', 3)
xline(918, 'Color', [0.5 0.85 0.5],'Linewidth', 3)
xline(629, 'Color', [0.5 0.85 0.5],'Linewidth', 3)
xline(493, 'Color', [0.5 0.85 0.5],'Linewidth', 3)
xline(834, 'Color', [0.5 0.85 0.5],'Linewidth', 3)
legend(cellstr(Regcoeff.i), 'Location', 'Eastoutside')

set(gcf,'renderer','Painters')
savename = ['C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/RSstrawb_regcoef_month' num2str(aopt)];
saveas(gcf,savename,'png')
saveas(gcf,savename,'epsc')
saveas(gcf,savename,'fig')




%% VALIDATION SCHEME  SUNNY SIDE -> SHADOW SIDE RAMAN
% -------------------------------------------------------------------------

% Merge 2021 and 2022 data blocks
X_als = saisir_rowmerge(X_prep_als_subsets{3},X_prep_als_subsets{6});
X_alssnv = saisir_rowmerge(X_prep_alssnv_subsets{3},X_prep_alssnv_subsets{6});
Y = saisir_rowmerge(Y_subsets{3},Y_subsets{6});

% Validate from identifier (sunny/shadow side)
ncomp = 20;
target = 1;
[FIG, Regcoeff, Performance, Residuals] = PLSR_validate_from_identifier3(X_alssnv,Y,target,5,ncomp,[1:3,9:14],'plotLimits',[1 30]); %#ok<ASGLU> % 5 Best choice Raman with test set. CHck later if this is different number when decided from crossval on cal set.

plot_regcoef(X_alssnv.v,squeeze(Regcoeff.d(1,:,:))',0.5, 'findPeaks','no','scaledPlot','yes')
%ylim([-0.25 0.15])
title(Regcoeff.i(1,:))
plot_regcoef(X_alssnv.v,squeeze(Regcoeff.d(2,:,:))',0.5, 'findPeaks','no','scaledPlot','yes')
title(Regcoeff.i(2,:))
%ylim([-0.25 0.15])

%% LEAVE-ONE-OUT CROSS VALIDATION, NIR
% -------------------------------------------------------------------------
% To show that Raman and NIR have similar results, so that it is a good
% case to compare the dependence on number of samples on.

% Merge 2021 and 2022 data blocks
X_snv = saisir_rowmerge(X_prep_snv_subsets{9},X_prep_snv_subsets{12});
% X_deriv2 = saisir_rowmerge(X_prep_deriv2_subsets{9},X_prep_deriv2_subsets{12}); 
% X_alssnv = saisir_rowmerge(X_prep_alssnv_subsets{9},X_prep_alssnv_subsets{12});% 
Y = saisir_rowmerge(Y_subsets{9},Y_subsets{12});

idpos = [1:3,9:14];
aopt = 20;
target = 1;
[Predictions,Regcoeff, Performance] = PLSR_CV_from_identifier(X_snv, Y, target, idpos, aopt,'PLotLimits',[1 20]);

figs = findobj('type','Figure');
FIG1 = figs(2);
FIG2 = figs(1);

saveas(FIG1,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/NIRstrawb_predvstarget_LooCV_rmseVSncomp','png')
saveas(FIG1,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/NIRstrawb_predvstarget_LooCV_rmseVSncomp','epsc')
saveas(FIG1,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/NIRstrawb_predvstarget_LooCV_rmseVSncomp','fig')
saveas(FIG2,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/NIRstrawb_predvstarget_LooCV_predvstarget','png')
saveas(FIG2,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/NIRstrawb_predvstarget_LooCV_predvstarget','epsc')
saveas(FIG2,'C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/NIRstrawb_predvstarget_LooCV_predvstarget','fig')



plot_regcoef(X_snv.v,squeeze(Regcoeff.d(:,:,12)),0.5, 'findPeaks','no')

%% VALIDATION 2021 -> 2022 MikroNIR // OR BY MONTH
% -------------------------------------------------------------------------

% Merge 2021 and 2022 data blocks
% X_deriv2 = saisir_rowmerge(X_prep_deriv2_subsets{9},X_prep_deriv2_subsets{12});
% X_deriv1 = saisir_rowmerge(X_prep_deriv1_subsets{9},X_prep_deriv1_subsets{12});;
% X_deriv2snv = saisir_rowmerge(X_prep_deriv2snv_subsets{9},X_prep_deriv2snv_subsets{12});
% X_deriv1snv = saisir_rowmerge(X_prep_deriv1snv_subsets{9},X_prep_deriv1snv_subsets{12});
X_snv = saisir_rowmerge(X_prep_snv_subsets{9},X_prep_snv_subsets{12});
Y = saisir_rowmerge(Y_subsets{9},Y_subsets{12});

% % Look only a side A/B, one at a time
% side = 'A'; % 'B'
% X_snv_A = select_from_identifier(X_snv,5,side); 
% Y_A = select_from_identifier(Y,5,side);

% Have to decide num components from calibration set first... 

% Validate from identifier (year)
aopt = 5; 
[FIG, Regcoeff, Performance] = PLSR_validate_from_identifier2(X_snv,Y,1, ...
                                9:14,aopt,[1:3,9:14], 'aoptAlg','WM',...
                                'aoptAlgMetric','RMSEP','aoptMode','calCV',...
                                'plotLimits', [1 20]); %#ok<ASGLU> %      %9:12(14) for year only % 5-6 Best choice Raman with test set. CHck later if this is different number when decided from crossval on cal set.
lgd = findobj('type','legend');
delete(lgd)
name4save = ['C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/NIRstrawb_predvstarget_month_amax',num2str(aopt)];
set(gcf,'renderer','Painters')
saveas(gcf,name4save,'png')
saveas(gcf,name4save,'epsc')
saveas(gcf,name4save,'fig')

% Plot performance summary
[FIG3, ah] = plot_allcombos_performance(Performance);
set(ah{1}, 'YLim',[0, 1.7])
set(ah{2}, 'YLim',[-2 2])
set(ah{3}, 'YLim',[0, 1.5])

name4save = ['C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/NIRstrawb_errors_month_amax',num2str(aopt)];
set(gcf,'renderer','Painters')
saveas(gcf,name4save,'png')
saveas(gcf,name4save,'epsc')
saveas(gcf,name4save,'fig')


plot_regcoef(Regcoeff.v, Regcoeff.d, 0.95,'plotColors',plot_colors,'findpeaks','no')
xline(1200,'LineWidth',30,'Color',[0.5 0.85, 0.5])
xline(1437, 'LineWidth',30,'Color',[0.5 0.85, 0.5])
legend(cellstr(Regcoeff.i),'Location','Eastoutside')
xlabel('Wavelength (nm)','FontSize',18)
%set(gca, 'Children', flipud(get(gca, 'Children')) )
name4save = ['C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/NIRstrawb_regcoeff_month_amax',num2str(aopt)];
set(gcf,'renderer','Painters')
saveas(gcf,name4save,'png')
saveas(gcf,name4save,'epsc')
saveas(gcf,name4save,'fig')





%% VALIDATION SCHEME  SUNNY SIDE -> SHADOW SIDE MikroNIR
% -------------------------------------------------------------------------

% Merge 2021 and 2022 data blocks
X_snv = saisir_rowmerge(X_prep_snv_subsets{9},X_prep_snv_subsets{12});
Y = saisir_rowmerge(Y_subsets{9},Y_subsets{12});

% Validate from identifier (sunny/shadow side)
aopt = 5;
[FIG, Regcoeff, Performance, Residuals] = PLSR_validate_from_identifier3(X_snv,Y,1,5,aopt,[1:3,9:14],'plotLimits',[1 30]); %#ok<ASGLU> % 5 Best choice Raman with test set. CHck later if this is different number when decided from crossval on cal set.

plot_regcoef(X_snv.v,squeeze(Regcoeff.d(1,:,:))',0.5, 'findPeaks','no')
title(Regcoeff.i(1,:))

plot_regcoef(X_snv.v,squeeze(Regcoeff.d(2,:,:))',0.5, 'findPeaks','no')
title(Regcoeff.i(2,:))


%% VALIDATION AS A FUNCTION OF NUMBER OF SAMPLES: RAMAN
% -------------------------------------------------------------------------

% Merge 2021 and 2022 data blocks
%X_als = saisir_rowmerge(X_prep_als_subsets{3},X_prep_als_subsets{6}); 
X_alssnv = saisir_rowmerge(X_prep_alssnv_subsets{3},X_prep_alssnv_subsets{6});
%X_alsnormpeak = saisir_rowmerge(X_prep_alsnormpeak_subsets{3},X_prep_alsnormpeak_subsets{6}); 
Y = saisir_rowmerge(Y_subsets{3},Y_subsets{6});

% % Look only a side A/B, one at a time
% side = 'A'; % 'B'
% X_alssnv_A = select_from_identifier(X_alssnv,5,side); 
% Y_A = select_from_identifier(Y,5,side); %#ok<NASGU>

rng default
ncompmax = 15;
[Result_summary_Raman_testrmsepcor, Result_Perdraw_Raman_testrmsepcor] = PLSR_montecarlo_nsamplesval(X_alssnv,Y,1,1000, [3:3:130], 9:12, ncompmax,'aoptAlgMetric','RMSEP_corr'); % ,'aoptAlgMetric','RMSEP_corr'
[Result_summary_Raman_testrmsep,Result_Perdraw_Raman_testrmsep] = PLSR_montecarlo_nsamplesval(X_alssnv,Y,1,1000, [3:3:130], 9:12, ncompmax,'aoptAlgMetric','RMSEP'); % ,'aoptAlgMetric','RMSEP_corr'
[Result_summary_RamanCV,Result_Perdraw_RamanCV] = PLSR_montecarlo_nsamplesvalCV(X_alssnv,Y,1,1000, [3:3:130], 9:12,[1:3,9:14],ncompmax); %


% Plot mean RMSEPs as a function of number of samples
Result_summary_Raman = Result_summary_RamanCV;

figure
hold on
p1 = plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(2,1,:)),'b-','LineWidth',3);
p2 = plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(6,1,:)),'b--');
p3 = plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(7,1,:)),'b--');
xlabel('Number of samples', 'FontSize',18)
ylabel('Median RMSEP', 'FontSize',18)
legend([p1(1) p2(1) p3(1)],'Median','5th percentile','95th percentile')
box off

figure
hold on
plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(2,8,:)),'b-','LineWidth',3)
plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(6,8,:)),'b--')
plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(7,8,:)),'b--')
xlabel('Number of samples', 'FontSize',18)
ylabel('Median RMSEP_{corr}', 'FontSize',18)
legend('Median','5th percentile','95th percentile')
box off

figure
hold on
plot(Result_summary_Raman.nsamples, squeeze(Result_summary_Raman.DescrStatsMetrics.d(2,6,:)),'b-','LineWidth',3)
plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(6,6,:)),'b--')
plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(7,6,:)),'b--')
xlabel('Number of samples', 'FontSize',18)
ylabel('Median LV', 'FontSize',18)
legend('Median','5th percentile','95th percentile')
box off

figure
hold on 
plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(2,7,:)),'b-','LineWidth',3)
plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(6,7,:)),'b--')
plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(7,7,:)),'b--')
xlabel('Number of samples', 'FontSize',18)
ylabel('Median R2P_{corr}', 'FontSize',18)
legend('Median','5th percentile','95th percentile')
box off

figure
hold on
plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(2,2,:)),'b-','LineWidth',3)
plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(6,2,:)),'b--')
plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(7,2,:)),'b--')
xlabel('Number of samples', 'FontSize',18)
ylabel('Median R2P', 'FontSize',18)
legend('Median','5th percentile','95th percentile')
box off

figure
hold on
plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(2,3,:)),'b-','LineWidth',3)
plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(6,3,:)),'b--')
plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(7,3,:)),'b--')
xlabel('Number of samples', 'FontSize',18)
ylabel('Median BIAS', 'FontSize',18)
yline(0,'--')
legend('Median','5th percentile','95th percentile')
box off

figure
hold on
plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(2,4,:)),'b-','LineWidth',3)
plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(6,4,:)),'b--')
plot(Result_summary_Raman.nsamples,  squeeze(Result_summary_Raman.DescrStatsMetrics.d(7,4,:)),'b--')
xlabel('Number of samples', 'FontSize',18)
ylabel('Median SLOPE', 'FontSize',18)
yline(1,'--')
legend('Median','5th percentile','95th percentile')
box off


% Same plots with ntot = 100 only (REDUCED PLOT)


%% VALIDATION AS A FUNCTION OF NUMBER OF SAMPLES: MikroNIR
% -------------------------------------------------------------------------

% Merge 2021 and 2022 data blocks
X_snv = saisir_rowmerge(X_prep_snv_subsets{9},X_prep_snv_subsets{12}); 
Y = saisir_rowmerge(Y_subsets{9},Y_subsets{12});

% % Look only a side A/B, one at a time
% side = 'A'; % 'B'
% X_alssnv_A = select_from_identifier(X_snv,5,side); 
% Y_A = select_from_identifier(Y,5,side); %#ok<NASGU>

ncompmax = 15;
[Result_summary_NIR_testrmsepcor, Result_Perdraw_NIR_testrmsepcor] = PLSR_montecarlo_nsamplesval(X_snv,Y,1,1000, [3:3:130], 9:12, ncompmax,'aoptAlgMetric', 'RMSEP_corr'); % ,'aoptAlgMetric','RMSEP_corr'
[Result_summary_NIR_testrmsep, Result_Perdraw_NIR_testrmsep]  = PLSR_montecarlo_nsamplesval(X_snv,Y,1,1000, [3:3:130], 9:12, ncompmax,'aoptAlgMetric', 'RMSEP'); % ,'aoptAlgMetric','RMSEP_corr'
[Result_summary_NIRCV, Result_Perdraw_NIR_rmsecv]  = PLSR_montecarlo_nsamplesvalCV(X_snv,Y,1,1000, [3:3:130], 9:12,[1:3,9:14],ncompmax); %

% Plot mean RMSEPs as a function of number of samples
Result_summary_NIR = Result_summary_NIRCV;

figure
hold on
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(2,1,:)),'b-','LineWidth',3)
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(6,1,:)),'b--')
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(7,1,:)),'b--')
xlabel('Number of samples', 'FontSize',18)
ylabel('Median RMSEP', 'FontSize',18)
legend('Median','5th percentile','95th percentile')
box off

figure
hold on
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(2,8,:)),'b-','LineWidth',3)
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(6,8,:)),'b--')
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(7,8,:)),'b--')
xlabel('Number of samples', 'FontSize',18)
ylabel('RMSEP_{corr}', 'FontSize',18)
legend('Median','5th percentile','95th percentile')
box off

figure
hold on
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(2,6,:)),'b-','LineWidth',3)
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(6,6,:)),'b--')
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(7,6,:)),'b--')
xlabel('Number of samples', 'FontSize',18)
ylabel('Median LV', 'FontSize',18)
legend('Median','5th percentile','95th percentile')
box off

figure
hold on
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(2,7,:)),'b-','LineWidth',3)
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(6,7,:)),'b--')
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(7,7,:)),'b--')
xlabel('Number of samples', 'FontSize',18)
ylabel('Median R2P_{corr}', 'FontSize',18)
legend('Median','5th percentile','95th percentile')
box off

figure
hold on
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(2,2,:)),'b-','LineWidth',3)
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(6,2,:)),'b--')
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(7,2,:)),'b--')
xlabel('Number of samples', 'FontSize',18)
ylabel('Median R2P', 'FontSize',18)
legend('Median','5th percentile','95th percentile')
box off

figure 
hold on
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(2,3,:)),'b-','LineWidth',3)
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(6,3,:)),'b--')
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(7,3,:)),'b--')
xlabel('Number of samples', 'FontSize',18)
ylabel('Median BIAS', 'FontSize',18)
yline(0,'--')
legend('Median','5th percentile','95th percentile')
box off

figure
hold on
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(2,4,:)),'b-','LineWidth',3)
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(6,4,:)),'b--')
plot(Result_summary_NIR.nsamples,  squeeze(Result_summary_NIR.DescrStatsMetrics.d(7,4,:)),'b--')
xlabel('Number of samples', 'FontSize',18)
ylabel('Median SLOPE', 'FontSize',18)
yline(1,'--')
legend('Median','5th percentile','95th percentile')
box off

%% COMPARE RAMAN AND NIR #SAMPLES
% -------------------------------------------------------------------------

% Aopt from test set:
for i = 1:size(Result_summary_Raman.avgPerformanceMetric.d,2)
   
    figure
    hold on
    p1 = plot(Result_summary_NIR_testrmsep.nsamples,  Result_summary_NIR_testrmsep.avgPerformanceMetric.d(:,i),'-o','Color',plot_colors(1,:));
    p2 = plot(Result_summary_Raman_testrmsep.nsamples,  Result_summary_Raman_testrmsep.avgPerformanceMetric.d(:,i),'-o','Color',plot_colors(2,:)) ;
    xlabel('Total number of samples', 'FontSize',18)
    ylabel(Result_summary_Raman.avgPerformanceMetric.v(i,:), 'FontSize',18)
    box off
    legend([p1(1) p2(1)], 'MikroNIR', 'Raman Kaiser')

end



% Aopt from CV vs test set and based on different metrics:
for i = 1:size(Result_summary_Raman.avgPerformanceMetric.d,2)
   
    figure
    hold on
    p1 = plot(Result_summary_NIRCV.nsamples,  Result_summary_NIR.avgPerformanceMetric.d(:,i),'-o','Color',plot_colors(1,:));
    p2 = plot(Result_summary_RamanCV.nsamples,  Result_summary_Raman.avgPerformanceMetric.d(:,i),'-*','Color',plot_colors(1,:)) ;

    p3 = plot(Result_summary_NIR_testrmsep.nsamples,  Result_summary_NIR_testrmsep.avgPerformanceMetric.d(:,i),'-o','Color',plot_colors(2,:));
    p4 = plot(Result_summary_Raman_testrmsep.nsamples,  Result_summary_Raman_testrmsep.avgPerformanceMetric.d(:,i),'-*','Color',plot_colors(2,:)) ;

    p5 = plot(Result_summary_NIR_testrmsepcor.nsamples,  Result_summary_NIR_testrmsepcor.avgPerformanceMetric.d(:,i),'-o','Color',plot_colors(3,:));
    p6 = plot(Result_summary_Raman_testrmsepcor.nsamples,  Result_summary_Raman_testrmsepcor.avgPerformanceMetric.d(:,i),'-*','Color',plot_colors(3,:)) ;

    xlabel('Total number of samples', 'FontSize',18)
    ylabel(Result_summary_RamanCV.avgPerformanceMetric.v(i,:), 'FontSize',18)
    box off
    legend([p1(1) p2(1) p3(1) p4(1) p5(1) p6(1)], 'MikroNIR CVaopt', 'Kaiser CVaopt', 'MikroNIR Test RMSEP aopt', 'Kaiser Test RMSEP aopt' , 'MikroNIR Test RMSEPcorr aopt', 'Kaiser Test RMSEPcorr aopt' )

end

% Aopt from CV vs test set and based on different metrics, standard
% deviation in metrics
for i = 1:size(Result_summary_Raman.avgPerformanceMetric.d,2)
   
    figure
    hold on
    p1 = plot(Result_summary_NIRCV.nsamples,  Result_summary_NIR.stdPerformanceMetric.d(:,i),'-o','Color',plot_colors(1,:));
    p2 = plot(Result_summary_RamanCV.nsamples,  Result_summary_Raman.stdPerformanceMetric.d(:,i),'-*','Color',plot_colors(1,:)) ;

    p3 = plot(Result_summary_NIR_testrmsep.nsamples,  Result_summary_NIR_testrmsep.stdPerformanceMetric.d(:,i),'-o','Color',plot_colors(2,:));
    p4 = plot(Result_summary_Raman_testrmsep.nsamples,  Result_summary_Raman_testrmsep.stdPerformanceMetric.d(:,i),'-*','Color',plot_colors(2,:)) ;

    p5 = plot(Result_summary_NIR_testrmsepcor.nsamples,  Result_summary_NIR_testrmsepcor.stdPerformanceMetric.d(:,i),'-o','Color',plot_colors(3,:));
    p6 = plot(Result_summary_Raman_testrmsepcor.nsamples,  Result_summary_Raman_testrmsepcor.stdPerformanceMetric.d(:,i),'-*','Color',plot_colors(3,:)) ;

    xlabel('Total number of samples', 'FontSize',18)
    ylabel(Result_summary_RamanCV.stdPerformanceMetric.v(i,:), 'FontSize',18)
    box off
    legend([p1(1) p2(1) p3(1) p4(1) p5(1) p6(1)], 'MikroNIR CVaopt', 'Kaiser CVaopt', 'MikroNIR Test RMSEP aopt', 'Kaiser Test RMSEP aopt' , 'MikroNIR Test RMSEPcorr aopt', 'Kaiser Test RMSEPcorr aopt' )

end


%% COMPARE RAMAN AND NIR BOOTSTRAP METHOD, FOR #SAMPLES = 102 
% (Performances stabilized for that large number. and still leaves room
% for variation in drawn samples for upper brix range samples)
% -------------------------------------------------------------------------


row_nsamp102 = 34;


% for i = 1:size(Result_summary_Raman.DescrStatsMetrics.d,2)
% 
%     figure
%     hold on
%     p1 = errorbar(1, Result_summary_NIR.DescrStatsMetrics.d(2,i,row_nsamp102), abs(Result_summary_NIR.DescrStatsMetrics.d(6,i,row_nsamp102)- Result_summary_NIR.DescrStatsMetrics.d(2,i,row_nsamp102)), abs(Result_summary_NIR.DescrStatsMetrics.d(7,i,row_nsamp102)-Result_summary_NIR.DescrStatsMetrics.d(2,i,row_nsamp102)),'o','Color',plot_colors(1,:));%,'EdgeColor', plot_colors(1,:));
%     p3 = errorbar(2, Result_summary_NIR_testrmsep.DescrStatsMetrics.d(2,i,row_nsamp102), abs(Result_summary_NIR_testrmsep.DescrStatsMetrics.d(6,i,row_nsamp102)- Result_summary_NIR_testrmsep.DescrStatsMetrics.d(2,i,row_nsamp102)), abs(Result_summary_NIR_testrmsep.DescrStatsMetrics.d(7,i,row_nsamp102)-Result_summary_NIR_testrmsep.DescrStatsMetrics.d(2,i,row_nsamp102)),'o','Color',plot_colors(2,:));%,'EdgeColor', plot_colors(1,:));
%     p5 = errorbar(3, Result_summary_NIR_testrmsepcor.DescrStatsMetrics.d(2,i,row_nsamp102), abs(Result_summary_NIR_testrmsepcor.DescrStatsMetrics.d(6,i,row_nsamp102)- Result_summary_NIR_testrmsepcor.DescrStatsMetrics.d(2,i,row_nsamp102)), abs(Result_summary_NIR_testrmsepcor.DescrStatsMetrics.d(7,i,row_nsamp102)-Result_summary_NIR_testrmsepcor.DescrStatsMetrics.d(2,i,row_nsamp102)),'o','Color',plot_colors(3,:));%,'EdgeColor', plot_colors(1,:));
% 
%     p2 = errorbar(4, Result_summary_Raman.DescrStatsMetrics.d(2,i,row_nsamp102), abs(Result_summary_Raman.DescrStatsMetrics.d(6,i,row_nsamp102)- Result_summary_Raman.DescrStatsMetrics.d(2,i,row_nsamp102)), abs(Result_summary_Raman.DescrStatsMetrics.d(7,i,row_nsamp102)-Result_summary_Raman.DescrStatsMetrics.d(2,i,row_nsamp102)),'*','Color',plot_colors(1,:));%,'EdgeColor', plot_colors(1,:));
%     p4 = errorbar(5, Result_summary_Raman_testrmsep.DescrStatsMetrics.d(2,i,row_nsamp102), abs(Result_summary_Raman_testrmsep.DescrStatsMetrics.d(6,i,row_nsamp102)- Result_summary_Raman_testrmsep.DescrStatsMetrics.d(2,i,row_nsamp102)), abs(Result_summary_Raman_testrmsep.DescrStatsMetrics.d(7,i,row_nsamp102)-Result_summary_Raman_testrmsep.DescrStatsMetrics.d(2,i,row_nsamp102)),'*','Color',plot_colors(2,:));%,'EdgeColor', plot_colors(1,:));
%     p6 = errorbar(6, Result_summary_Raman_testrmsepcor.DescrStatsMetrics.d(2,i,row_nsamp102), abs(Result_summary_Raman_testrmsepcor.DescrStatsMetrics.d(6,i,row_nsamp102)- Result_summary_Raman_testrmsepcor.DescrStatsMetrics.d(2,i,row_nsamp102)), abs(Result_summary_Raman_testrmsepcor.DescrStatsMetrics.d(7,i,row_nsamp102)-Result_summary_Raman_testrmsepcor.DescrStatsMetrics.d(2,i,row_nsamp102)),'*','Color',plot_colors(3,:));%,'EdgeColor', plot_colors(1,:));
% 
% 
%     %xlabel('Method', 'FontSize',18)
%     ylabel(Result_summary_RamanCV.DescrStatsMetrics.v(i,:), 'FontSize',18)
%     box off
%     %legend([p1(1) p2(1) p3(1) p4(1) p5(1) p6(1)], 'MikroNIR CVaopt', 'Kaiser CVaopt', 'MikroNIR Test RMSEP aopt', 'Kaiser Test RMSEP aopt' , 'MikroNIR Test RMSEPcorr aopt', 'Kaiser Test RMSEPcorr aopt' )
%     set(gca,'XTickLabel',{'MikroNIR_{CVaopt}',  'MikroNIR_{RMSEPaopt}', 'MikroNIR_{RMSEPCaopt}','Kaiser_{CVaopt}', 'Kaiser_{RMSEPaopt}' , 'Kaiser_{RMSEPCaopt}'})
%     xtickangle(90)
% 
% end


% Violin plots ------------------------------------------------------------
errortype = 3; %8,6,4,3
AOPTMAX = input('Enter the maxium component restriction: ');

ViolinBootTestRMSEPcorr = Result_Perdraw_Raman_testrmsepcor.d(34,:,errortype); % TOT SAMPLES 102, ALL DRAWS, RMSEP (AOPT)
ViolinBootTestRMSEP = Result_Perdraw_Raman_testrmsep.d(34,:,errortype);
ViolinBootRMSECV = Result_Perdraw_RamanCV.d(34,:,errortype);
NIRViolinBootTestRMSEPcorr = Result_Perdraw_NIR_testrmsepcor.d(34,:,errortype);
NIRViolinBootTestRMSEP = Result_Perdraw_NIR_testrmsep.d(34,:,errortype);
NIRViolinBootRMSECV = Result_Perdraw_NIR_rmsecv.d(34,:,errortype);

BootStrapModes = {'Raman_{RMSEPcorr}';'Raman_{RMSEP}';'Raman_{RMSECV}';...
                  'NIR_{RMSEPcorr}';'NIR_{RMSEP}';'NIR_{RMSECV}'};

figure
vs = violinplot([ViolinBootTestRMSEPcorr' ViolinBootTestRMSEP' ViolinBootRMSECV' NIRViolinBootTestRMSEPcorr' NIRViolinBootTestRMSEP' NIRViolinBootRMSECV'],...
                BootStrapModes,'MarkerSize', 1,'MedianMarkerSize', 45, 'QuartileStyle','shadow');
ylabel(Result_Perdraw_NIR_rmsecv.k(errortype,:))
%ylim([0.5 1.5])

% ylim([0 1.2])
% yline(1)

% ylim([-4 5])
% yline(0)


box off 
grid on
set(gcf,'renderer','Painters')

etype = Result_Perdraw_NIR_rmsecv.k(errortype,:);
saveas(gcf,['C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/ViolinplotBootstrap_Amax',num2str(AOPTMAX),'_',etype(1:end-9)],'png')
saveas(gcf,['C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/ViolinplotBootstrap_Amax',num2str(AOPTMAX),'_',etype(1:end-9)],'epsc')
saveas(gcf,['C:\Users\Tiril Aurora\Documents\01 Nofima\MATLAB\NIR vs Raman robustness paper\Illustrations/ViolinplotBootstrap_Amax',num2str(AOPTMAX),'_',etype(1:end-9)],'fig')


% % Aopt from CV vs test set and based on different metrics:
% for i = 1:size(Result_summary_Raman.avgPerformanceMetric.d,2)
% 
%     figure
%     hold on
%     p1 = bar({'MikroNIR CVaopt'},  Result_summary_NIR.avgPerformanceMetric.d(row_nsamp102,i),'FaceColor',plot_colors(1,:),'EdgeColor', plot_colors(1,:));
%     p3 = bar({'MikroNIR Test RMSEP aopt'},  Result_summary_NIR_testrmsep.avgPerformanceMetric.d(row_nsamp102,i),'FaceColor',plot_colors(2,:),'EdgeColor', plot_colors(2,:));
%     p5 = bar({'MikroNIR Test RMSEPcorr aopt'},  Result_summary_NIR_testrmsepcor.avgPerformanceMetric.d(row_nsamp102,i),'FaceColor',plot_colors(3,:),'EdgeColor', plot_colors(3,:));
% 
%     p2 = bar({'Raman Kaiser CVaopt'},  Result_summary_Raman.avgPerformanceMetric.d(row_nsamp102,i),'FaceColor',plot_colors(1,:),'EdgeColor', plot_colors(1,:)) ;
%     p4 = bar({'Raman Kaiser Test RMSEP aopt'},  Result_summary_Raman_testrmsep.avgPerformanceMetric.d(row_nsamp102,i),'FaceColor',plot_colors(2,:),'EdgeColor', plot_colors(2,:)) ;
%     p6 = bar({'Raman Kaiser Test RMSEPcorr aopt'},  Result_summary_Raman_testrmsepcor.avgPerformanceMetric.d(row_nsamp102,i),'FaceColor',plot_colors(3,:),'EdgeColor', plot_colors(3,:)) ;
% 
%     %xlabel('Method', 'FontSize',18)
%     ylabel(Result_summary_RamanCV.avgPerformanceMetric.v(i,:), 'FontSize',18)
%     box off
%     %legend([p1(1) p2(1) p3(1) p4(1) p5(1) p6(1)], 'MikroNIR CVaopt', 'Kaiser CVaopt', 'MikroNIR Test RMSEP aopt', 'Kaiser Test RMSEP aopt' , 'MikroNIR Test RMSEPcorr aopt', 'Kaiser Test RMSEPcorr aopt' )
% 
% end


% Aopt from CV vs test set and based on different metrics, with error bar (std):
% NOTE: THE STD COULD BE INTERCHANGED WITH THE PERCENTILE VALUES (2.5% -
% 97,5%) TO OBTAIN 95% CONFIDENCE INTERVAL. NOTE THAT THIS IS APPROXIMATELY
% THE SAME AS +/- 2*STD
% for i = 1:size(Result_summary_Raman.avgPerformanceMetric.d,2)
% 
%     figure
%     hold on
%     p1 = errorbar(1, Result_summary_NIR.avgPerformanceMetric.d(row_nsamp102,i), Result_summary_NIR.stdPerformanceMetric.d(row_nsamp102,i),'o','Color',plot_colors(1,:));%,'EdgeColor', plot_colors(1,:));
%     p3 = errorbar(2,  Result_summary_NIR_testrmsep.avgPerformanceMetric.d(row_nsamp102,i), Result_summary_NIR_testrmsep.stdPerformanceMetric.d(row_nsamp102,i),'o','Color',plot_colors(2,:));%,'EdgeColor', plot_colors(2,:));
%     p5 = errorbar(3,  Result_summary_NIR_testrmsepcor.avgPerformanceMetric.d(row_nsamp102,i), Result_summary_NIR_testrmsepcor.stdPerformanceMetric.d(row_nsamp102,i),'o','Color',plot_colors(3,:));%,'EdgeColor', plot_colors(3,:));
% 
%     p2 = errorbar(4,  Result_summary_Raman.avgPerformanceMetric.d(row_nsamp102,i), Result_summary_Raman.stdPerformanceMetric.d(row_nsamp102,i),'*','Color',plot_colors(1,:));%,'EdgeColor', plot_colors(1,:)) ;
%     p4 = errorbar(5,  Result_summary_Raman_testrmsep.avgPerformanceMetric.d(row_nsamp102,i),  Result_summary_Raman_testrmsep.stdPerformanceMetric.d(row_nsamp102,i),'*','Color',plot_colors(2,:));%,'EdgeColor', plot_colors(2,:)) ;
%     p6 = errorbar(6,  Result_summary_Raman_testrmsepcor.avgPerformanceMetric.d(row_nsamp102,i), Result_summary_Raman_testrmsepcor.stdPerformanceMetric.d(row_nsamp102,i),'*','Color',plot_colors(3,:));%,'EdgeColor', plot_colors(3,:)) ;
% 
%     %xlabel('Method', 'FontSize',18)
%     ylabel(Result_summary_RamanCV.avgPerformanceMetric.v(i,:), 'FontSize',18)
%     box off
%     %legend([p1(1) p2(1) p3(1) p4(1) p5(1) p6(1)], 'MikroNIR CVaopt', 'Kaiser CVaopt', 'MikroNIR Test RMSEP aopt', 'Kaiser Test RMSEP aopt' , 'MikroNIR Test RMSEPcorr aopt', 'Kaiser Test RMSEPcorr aopt' )
%     set(gca,'XTickLabel',{'MikroNIR_{CVaopt}',  'MikroNIR_{RMSEPaopt}', 'MikroNIR_{RMSEPCaopt}','Kaiser_{CVaopt}', 'Kaiser_{RMSEPaopt}' , 'Kaiser_{RMSEPCaopt}'})
%     xtickangle(90)
% 
% end
% 


