function [peaks,means,widths,fitsvalue,fitresult, gof] = fitNGaussains(N,x,y,peakcut,plotting)

if nargin==4
    plotting=0;
end
%x(y<=0)=[];
%y(y<=0)=[];
maxy=max(y);
y=y/maxy;
%y=y-median(y); % median baseline subtracte

minpeakheight=max(y)*peakcut;
[peakEstimate, peakXInd]=findpeaks(y, 'NPEAKS',N,'MINPEAKDISTANCE',3,'MINPEAKHEIGHT',minpeakheight); %median(y)*1.5

peakPosEstimate=x(peakXInd);


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
%ft = fittype( 'a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + a3*exp(-((x-b3)/c3)^2) + a4*exp(-((x-b4)/c4)^2) + a5*exp(-((x-b5)/c5)^2) + a6*exp(-((x-b6)/c6)^2) + a7*exp(-((x-b7)/c7)^2) + a8*exp(-((x-b8)/c8)^2) + a9*exp(-((x-b9)/c9)^2) + a10*exp(-((x-b10)/c10)^2) + a11*exp(-((x-b11)/c11)^2) + a12*exp(-((x-b12)/c12)^2) + a13*exp(-((x-b13)/c13)^2) + a14*exp(-((x-b14)/c14)^2) + a15*exp(-((x-b15)/c15)^2) + a16*exp(-((x-b16)/c16)^2) + a17*exp(-((x-b17)/c17)^2) + a18*exp(-((x-b18)/c18)^2) + aN*exp(-((x-bN)/cN)^2) + d', 'independent', 'x', 'dependent', 'y' );
ft = fittype([makeNGaussainEqaution(N)]);  % '+ p2*x^2 + p3*x'
opts = fitoptions( ft );
opts.Display = 'Off';

opts.Lower = [zeros(1,N) peakPosEstimate'-2 ones(1,N)*1 0 ]; %-inf -inf 

opts.StartPoint = [peakEstimate' peakPosEstimate' ones(1,N)*2 0.3]; %-2e-5 0.003

opts.Upper = [peakEstimate'*1.2 peakPosEstimate'+2 ones(1,N)*10 1]; % inf inf


% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

fitsvalue=coeffvalues(fitresult);
peaks=fitsvalue(1:N);
means=fitsvalue(N+1:N*2);
widths=fitsvalue(N*2+1:end-1);

if gof.rsquare<0.95
    warning('N gauss fit is not that great, inspect manually')
    plotting=1;
end

if plotting
    figure(plotting);
    subplot(2,1,1)
    plot(x,y,peakPosEstimate,peakEstimate,'o')
    %% Plot fit with data.
    subplot(2,1,2)
    h = plot( fitresult, xData, yData,'x');
    legend( h, 'y vs. x', 'untitled fit 1', 'Location', 'NorthEast' );
    % Label axes
    xlabel( 'x' );
    ylabel( 'y' );
    grid on
    hold on 
    line([0 max(x)], [minpeakheight minpeakheight]);
    fitresult
    gof
    pause;
end

