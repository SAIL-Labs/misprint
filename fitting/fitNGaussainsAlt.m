function [peaks,means,widths,xfitted] = fitNGaussainsAlt(N,x,y,peakcut,plotting)

if nargin==4
    plotting=0;
end
x(y<=0)=[];
y(y<=0)=[];
maxy=max(y);
y=y/maxy;
%y=y-median(y); % median baseline subtracte

minpeakheight=max(y)*peakcut;
[peakEstimate, peakXInd]=findpeaks(y, 'NPEAKS',N,'MINPEAKDISTANCE',3,'MINPEAKHEIGHT',minpeakheight); %median(y)*1.5

peakPosEstimate=x(peakXInd);


%% fit
options = optimset('Display','off',...
    'Algorithm','levenberg-marquardt',...
    'TolFun',1e-8,...
    'TolX',1e-8);

x0=[peakEstimate' peakPosEstimate' ones(1,N)*2 0.1];


fun = @(co,xData) sum(nGausFunc(co,xData,N),2);

[xfitted,resnorm] = lsqcurvefit(fun,x0,x,y,[],[],options);

%% export
peaks=xfitted(1:N)*maxy;
xfitted(1:N)=peaks;
means=xfitted(N+1:N*2);
widths=xfitted(N*2+1:end-1);
xfitted(end)=xfitted(end)*maxy;

%% checking
if resnorm > 0.4
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

