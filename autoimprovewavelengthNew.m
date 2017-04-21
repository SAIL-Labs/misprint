function autoimprovewavelengthNew(filename,varargin)
import chrislib.fitting.*
%filename='thar-0002-1D-spectra.fits';

spectra=fitsread(filename);
numOrders=size(spectra,3);
numFibres=size(spectra,1);

if isempty(varargin)
    fibreToUse=10;
else
    fibreToUse=varargin{1};
end

for o=1:numOrders
    spectra(:,:,o)=bsxfun(@rdivide,spectra(:,:,o),max(spectra(:,:,o)')');
end

orders=1:numOrders;
for o=orders
    matpayload=load([stripextension(filename) '-ref-points-order' num2str(o) '.mat']);
    pin(fibreToUse,:,o)=matpayload.p;
    xDatainit(o)={matpayload.xData};
    xRef(o)={matpayload.xRef};
    wavefit(fibreToUse,:,o)=polyval(pin(fibreToUse,:,o),1:size(spectra,2));
    
    IN1=spectra(fibreToUse,:,o)/max(spectra(fibreToUse,:,o));
    for f=[1:size(spectra,1)]
        IN2=spectra(f,:,o)/max(spectra(f,:,o));
        [lags(f,:, o),C(f,:, o),Info(f)]=xcorr_fft(IN1',IN2');
        shift(f, o)=Info(f).BestShift;
    end
end

%%
x=1:2504;
win=15;
for o=orders
    xDataFit=[];
    xData=[];
    xRefCur=xRef{o}';
    for f=1:numFibres
        xData(:,f)=xDatainit{o}-shift(f,o);
        for x=1:size(xData,1)
            xprofile=max([1 round(xData(x,f))-win]):min([size(spectra,2) round(xData(x,f))+win]);
            if length(xprofile)<20
                continue
            end
            xzoomprofile=spectra(f,xprofile,o);
            xzoomprofile=xzoomprofile/max(xzoomprofile);
            [cf_, b, c, d, a, gof] = fit_gauss(xprofile',xzoomprofile');
            
            if gof.rsquare<0.75
%                 xx=linspace(min(xprofile),max(xprofile),300);
%                 plot(xprofile,xzoomprofile,b,1,'rx',xx,feval(cf_,xx))
%                 title(gof.rsquare)
%                 hold on           
            else
                xDataFit(x,f)=b;
            end
        end
    end
    %hold off
    sum(~any(xDataFit==0,2))
    xRefCur(any(xDataFit==0,2))=[];
    xDataFit(any(xDataFit==0,2),:)=[];
    
    %%
    for f=1:numFibres
        [p(f,:,o),S(f,:,o),mu(f,:,o)] = polyfit(xDataFit(:,f),xRefCur,3);
    end
    save([stripextension(filename) '-autoFittedWave.mat'],'p','S','mu','shift')
end