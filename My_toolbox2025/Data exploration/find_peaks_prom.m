function [xindices,xvals] = find_peaks_prom(X,wn,row,threshold,gap,min_max)
%find_peaks   - find peaks greater than a threshold value 
% inside a window of size (data points) defined by gap. gap is a odd number
% return the VECTOR of positions (name of variable converted into numbers) 
% function vect=find_peaks(saisir,row,threshold,gap,(min_max))
% if min_max is larger than 0, the minimum are also identified. 
%
% wn - char array (as from saisir structure)
%
% Suggested default values for baseline corr spectra: vect = find_peaks(X,20,200,10,0)
% Adapted from EMSC toolox (Biospec group)
%
% Note: does not work properly if peaks are in oppsite direction 
%       (as could be the case in principal comps and reg. coefs etc.) It
%       then identifies the "minimums" in negative weights. USE: find_peaks_pca 

xvals=[];
[n,p]=size(X);
%peak.d=zeros(n,8);
if(nargin<5)
    min_max=0;
end
nel=0;
for i=row
   %disp(i);
   %courbe(X,i);
   x=X(i,:);
   aux=floor((gap+1)/2);
   aux1=aux-1;
   for j=aux:p-aux
      if(abs(x(j))>threshold)
         [a,b]=sort(-x(j-aux1:j+aux1));
         if (b(1)==aux)
            posx=str2num(wn(j,:));
            nel=nel+1;
            xvals(nel)=posx;
            xindices(nel) = j;
            %line([posx posx],[0 x(j)]);
         	%text(posx,x(j),num2str(round(str2num(wn(j,:)),1)),'FontSize',8);   
         end   
      end  
    end   
 end
vect1=xvals;


if(min_max>0)
clear xvals;
nel=0;
for i=row
   %disp(i);
   %courbe(saisir,i);
   x=X(i,:);
   aux=(gap+1)/2;
   aux1=aux-1;
   for j=aux:p-aux
      if(abs(x(j))>threshold)
      [a,b]=sort(-x(j-aux1:j+aux1));
      %   [b(end) aux]
         if (b(end)==aux)
            'ici'
            posx=str2num(wn(j,:));
            nel=nel+1;
            xvals(nel)=posx;
            %line([posx posx],[0 x(j)]);
            %text(posx,x(j),num2str(round(str2num(wn(j,:)),1)),'FontSize',8);   
         end   
      end  
    end   
 end
 xvals=[vect1 xvals];
end