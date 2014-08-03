%auto improove wavelength solution
%function [estimates start_point]=autoimprovewavelength
filename='43-1D-spectra';
spectra=fitsread([filename '.fits']);
for or=1:13
    spectra(:,:,or)=bsxfun(@rdivide,spectra(:,:,or),max(spectra(:,:,or)')');
end


for o=1:size(spectra,3)
    matpayload=load([filename '-ref-points-order' num2str(5) '.mat']);
    p(10,:,o)=matpayload.p;
    xDatainit(o)={matpayload.xData};
    xRef(o)={matpayload.xRef};
    wavefit(10,:,o)=polyval(p(10,:,o),1:size(spectra,2));
end


for order=1:size(spectra,3)
    IN1=spectra(10,:,order)/max(spectra(10,:,order));
    for fibre=[1:9 11:size(spectra,1)]
        IN2=spectra(fibre,:,order)/max(spectra(fibre,:,order));
        [lags(fibre,:, order),C(fibre,:, order),Info(fibre)]=xcorr_fft(IN1',IN2');
        shift(fibre,:, order)=Info(fibre).BestShift;
    end
end




%% find centroids of xref

for o=1;
    for f=1:19;
        dataspectra=(spectra(f,:,o));
        
%         figure(2)
%         plot(xDatainit{o}-shift(f,:, o),0.25,'rx',1:2537,spectra(f,:,o))
        
        xfit=[];
        for i=1:length(xDatainit{o})
            %x=[round(xDatainit{o}(i)-shift(f,:, o))-10:round(xDatainit{o}(i)-shift(f,:, o))+10];
            x=max([1 round(xDatainit{o}(i)-shift(f,:, o))-10]):min([2537 round(xDatainit{o}(i)-shift(f,:, o))+10]);
            xzoomprofile=dataspectra(x);
            
            [~,ind]=max(xzoomprofile);
            ind=ind-11;
            x=max([1 round(xDatainit{o}(i)-shift(f,:, o))-10+ind]):min([2537 round(xDatainit{o}(i)-shift(f,:, o))+10+ind]);
            xzoomprofile=dataspectra(x);
            xzoomprofile=xzoomprofile-xzoomprofile(1);
            
            [cf_, b, c, d, a, gof] = fit_gauss(x',xzoomprofile');
            xfit(i)=b(a==max(a));
             
            if gof.rsquare<0.95
                figure(1);clf
                plot(cf_,x,xzoomprofile,'x')
                hold all
                plot(xfit(i),a(a==max(a)),'og')
                hold off
                gof
                pause;
            end
        end
        
       
        
        xDatafit(o,f)={xfit};
        p(f,:,o)=polyfit(xDatafit{o,f},xRef{o},3);
    end
    
end

%%
wavelinear=[];
for o=1;
    for f=1:19;
        wavefit(f,:,o)=polyval(p(f,:,o),1:2537);
    end
    wavelinear(:,o)=linspace(min(wavefit(10,:,o)),max(wavefit(10,:,o)),2537);
end

speclinear=[];
for o=1;
    for f=1:19;
        speclinear(f,:,o)=interp1(wavefit(f,:,o),spectra(f,:,o),wavelinear(:,o));
    end
end
%plot(wavefit(:,:,1)', spectra(:,:,1)')

plot(wavelinear(:,1),speclinear(:,:,1))

wavelinear=wavelinear(~isnan(sum(speclinear,1)),:);
speclinear=speclinear(:,~isnan(sum(speclinear,1)),:);

clear lags shift C Info
for order=1:size(speclinear,3)
    IN1=speclinear(10,:,order)/max(speclinear(10,:,order));
    for fibre=[1:9 11:size(speclinear,1)]
        IN2=speclinear(fibre,:,order)/max(speclinear(fibre,:,order));
        [lags(fibre,:, order),C(fibre,:, order),Info(fibre)]=xcorr_fft(IN1',IN2');
        shift(fibre,:, order)=Info(fibre).BestShift*mean(diff(wavelinear));
    end
end



















return
% %%
%
%
%
%
%
%
% %options.funcCount=1e5;
% options = psoptimset(@patternsearch);
% options.TolMesh=1e-13;
% options.TolX=1e-13;
% options.TolFun=1e-13
% options.MaxIter=4000
% o=1;
% for f=8
%
%
%     start_point=p(o,:);
%     start_point(4)=start_point(4)+start_point(3)*shift(f,:, o);
%
%     [THWave,THIntensity] = ThArSpec(min2(wavefit(10,:,1))*10-10,max2(wavefit(10,:,1))*10+10);
%     THIntensity=THIntensity/max(THIntensity);
%     THWave=THWave/10;
%
%     model = @(x) expfun(x,spectra(f,:,o),THIntensity,THWave);
%
%     LB=start_point-abs(start_point).*[0.01 0.01 0.01 0.001];
%     UB=start_point+abs(start_point).*[0.01 0.01 0.01 0.001];
%
%     estimates(f,:,o) = patternsearch(model, start_point,[],[],[],[],start_point-abs(start_point).*0.03,start_point+abs(start_point)*0.03,[],options);
%
%     %estimates(f,:,o) =  simulannealbnd(model, start_point,start_point-abs(start_point).*0.1,start_point+abs(start_point)*0.1);
%     figure(f);clf
%     plot(polyval(estimates(f,:,o),1:2537),spectra(f,:,o),'b',polyval(start_point,1:2537),spectra(f,:,o),'g--',THWave,THIntensity,'r')
% end
% %error(' ')
% return
%
% figure(20)
%
% for i=1:19
%     plot(polyval(estimates(i,:,1),1:2537),spectra(i,:,1))
%     hold all
% end
% hold off
%
% end
% function [sse] = expfun(params,spectra,THIntensity,THWave)
% wavestart=polyval(params,1:2537);
% interpspec=interp1(wavestart,spectra,THWave,'linear',0);
% %[c]=xcorr(THIntensity,interpspec,'coeff');
%
% %sse=1-max(c)
% ErrorVector = interpspec - THIntensity;
% sse = sum(ErrorVector .^ 2);
% end
