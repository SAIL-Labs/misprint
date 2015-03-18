function [peaks,means,widths,xfitted] = fitNGaussainsAlt(N,x,y,peakcut,plotting)
import chrislib.*
if nargin==4
    plotting=0;
end
x(y<=0)=[];
y(y<=0)=[];
maxy=max(y);
y=y/maxy;

y=y-mean(y(1:4));
%y=y-median(y); % median baseline subtracte

minpeakheight=max(y)*peakcut;
if N==1
    peakXInd=find(max(y)==y,1,'first');
    peakEstimate=y(peakXInd);
else
    [peakEstimate, peakXInd]=findpeaks(y, 'NPEAKS',N,'MINPEAKDISTANCE',2,'MINPEAKWIDTH',2,'MinPeakProminence',0.1); %median(y)*1.5
end
peakPosEstimate=x(peakXInd);


%% fit
options = optimset('Display','off',...
    'TolFun',1e-8,...
    'TolX',1e-8);
%    'Algorithm','levenberg-marquardt',...

x0=[peakEstimate' peakPosEstimate' ones(1,N) 0.1];
xlb=[peakEstimate'*0.8 peakPosEstimate'*0.99 ones(1,N)*0.8 0];
xub=[peakEstimate'*1.2 peakPosEstimate'*1.01 ones(1,N)*3 0.2];


fun = @(co,xData) sum(misprint.nGausFunc(co,xData,N),2);

[xfitted,resnorm] = lsqcurvefit(fun,x0,x,y,xlb,xub,options);

%% export
peaks=xfitted(1:N)*maxy;
xfitted(1:N)=peaks;
means=xfitted(N+1:N*2);
widths=xfitted(N*2+1:end-1);
xfitted(end)=xfitted(end)*maxy;

%% checking
if resnorm > 2
    warning('N gauss fit is not that great, inspect manually')
    plotting=1;
end

if plotting
    figure(plotting);
    subplot(2,1,1)
    plot(x,y,peakPosEstimate,peakEstimate,'o')
    %% Plot fit with data.
    subplot(2,1,2)
    h = plot(x,fun(xfitted,x), x, y*maxy,'x');
    legend( h, 'y vs. x', 'untitled fit 1', 'Location', 'NorthEast' );
    % Label axes
    xlabel( 'x' );
    ylabel( 'y' );
    grid on
    hold on
    line([min(x) max(x)], [minpeakheight minpeakheight]*maxy);
    disp(['resnorm: ' num2str(resnorm)])
    hold off
    pause;
end

