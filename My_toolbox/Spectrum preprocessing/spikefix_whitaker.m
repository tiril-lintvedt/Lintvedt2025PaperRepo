function [x_fixed,spikes] = spikefix_whitaker(x, d, w, threshold, plotcheck)
% ------------------- Spike detection for Raman spectra -------------------  
%     Detect spikes based on the modified Z score, and removes the detected 
%     spikes by interpolation. Whitaker and Hayes propose to make advantage 
%     of the high intensity and small width of spikes and therefore use the 
%     difference between consecutive spectrum points ∇x(i) = x(i)-x(i-1) 
%     to calculate the z-scores, where x(1), …, x(n) are the values of a
%     single Raman spectrum.
%
%
%   INPUT:
%                       y - Spectrum
%                       m - integer, defines that we use 2m + 1 points around  
%                           the spike to compute and average which replaces 
%                           spike
%                       d - integer, degree of derivative to apply
%                       w - integer, buffer for spike width
%               threshold - Modified Z value threshold for detection of spikes
%                    plot - 0 or 1, whether to make inspection plots
% 
%   OUTPUT:
%                x_fixed - Spectrum without spikes
% 
%   Reference: A Simple Algorithm for Despiking Raman Spectra, Whitaker & Hayes, 2018 
%              python code: https://towardsdatascience.com/removing-spikes-from-raman-spectra-8a9fdda0ac22
%
%   Modifications: Use 2nd derivative in stead of 1st derivative (e.g lipid
%                  Raman bands are very strong/steep and comparable to spikes in 
%                  1st diff)
% 
%   Reported bugs: 
%                   If spike is detected in end regions of spectrum, then the spike-
%                   removing algorithm is will try to access out-of-bounds array 
%                   indices/elements.
% -------------------------------------------------------------------------

x_fixed = x;
spikes = abs(modified_z_score(diff(x,d))) > threshold ; % detect spikes using modified Z score

% Spike buffer: how wide of a spike do we want to account for ? 

% sp = find(spikes);
% % Spike width to left
% wleft = 1;
% for ip = 1:wleft
%     spikes(sp - ip) = 1;
% end
% 
% % Spike width to right
% wright = 2;
% for ip = 1:wright
%     spikes(sp + ip) = 1;
% end

sp = find(spikes);
%w = 2;
for ip = 1:w
    spikes(sp - ip) = 1;
    spikes(sp + ip) = 1;
end


 for i = 1:length(spikes)
     if spikes(i) ~= 0            % if we have a spike at position i
 
        % Find first point on each side of the spike, which is not a spike.
        jright = 1; jleft = 1;
        while spikes(i-jleft) == 1  % Can be unfortunate to get spike in position 1. Indexing not possible below 1.            
            jleft = jleft +1;
        end
        
        while spikes(i+jright) == 1
           jright = jright+1; 
        end
        
        w = (i-jleft):(i+jright); % Use this interval for spike correction
        x_fixed = replace_by_line(x_fixed, w);  % replace spike  with fitted line   
        i = i + jright; % Skip points which was fitted in this iteration

%        w = (i-m):(i+1+m);        % select 2m + 1 points around the spike 
%        % Make sure that this interval contains non-spike points
%        j = 1;
%        while sum(spikes(w) == 0) == 0 
%            w = (i-m-j):(i+1+m+j);
%            j = j+1;
%        end
%        
%        w2 = w(spikes(w) == 0);   % from this interval, choose points which are not spikes
%        x_fixed(i) = mean(x(w2)); % replace spikes with mean of surrounding non-spike points
        
     end
    

 end
 
% Inspection plots:
if plotcheck
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/3 50 scrsz(3)/2 scrsz(4)-150])
    subplot(7,1,1); plot(x); ylabel('original')
    subplot(7,1,2); plot(spikes); ylabel('selected spikes')
    subplot(7,1,3); plot(x_fixed); ylabel('fixed')
    subplot(7,1,4); plot(diff(x)); ylabel('1st diff')
    subplot(7,1,5); plot(modified_z_score(diff(x,1))); ylabel('mod. z (diff 1)')
    subplot(7,1,6); plot(diff(x,2)); ylabel('2nd diff')
    subplot(7,1,7); plot(modified_z_score(diff(x,d))); ylabel('mod. z (diff used)') 
end
end
