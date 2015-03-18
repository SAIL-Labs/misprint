%auto improove wavelength solution
function autoimprovewavelength(varargin)
    import chrislib.fitting.*
    
    if ~isempty(varargin)
        filename=varargin{1};
    else
        filename='arc-001-reduced-1D-spectra';
    end
    
    spectra=(fitsread([filename '.fits']));
    
    if length(varargin)==2
        newspectra(1,:,:)=spectra';
        spectra=newspectra; % flip orders and fibers for 1 core case
    end
    
    numOrders=size(spectra,3);
    
    for or=1:numOrders
        spectra(:,:,or)=bsxfun(@rdivide,spectra(:,:,or),max(spectra(:,:,or)')');
    end
    
    spectra(isnan(spectra))=0;
    try
        load([filename '-autoFittedWave.mat '])
    catch err
        
        for o=1:numOrders
            matpayload=load([filename '-ref-points-order' num2str(o) '.mat']);
            p(1,:,o)=matpayload.p;
            xDatainit(o)={matpayload.xData};
            xRef(o)={matpayload.xRef};
            wavefit(1,:,o)=polyval(p(1,:,o),1:size(spectra,2));
        end
        
        
        for order=1:numOrders
            IN1=spectra(1,:,order)/max(spectra(1,:,order));
            for fibre=[1:size(spectra,1)]
                IN2=spectra(fibre,:,order)/max(spectra(fibre,:,order));
                [lags(fibre,:, order),C(fibre,:, order),Info(fibre)]=xcorr_fft(IN1',IN2');
                shift(fibre,:, order)=Info(fibre).BestShift;
            end
        end
        
        avgshift=median(shift,3);
        
        
        %% find centroids of xref
        
        win=6;
        for o=1:numOrders;
            for f=1:size(spectra,1);
                disp(['Order: ' num2str(o) '  Fiber: ' num2str(f)]);
                dataspectra=spectra(f,:,o);
                
                %         figure(2)
                %         plot(xDatainit{o}-shift(f,:, o),0.25,'rx',1:2537,spectra(f,:,o))
                
                xfit=[];
                usemask=true(size(xDatainit{o}));
                for i=1:length(xDatainit{o})
                    try
                        %x=[round(xDatainit{o}(i)-shift(f,:, o))-10:round(xDatainit{o}(i)-shift(f,:, o))+10];
                        x=max([1 round(xDatainit{o}(i)-avgshift(f))-win]):min([size(spectra,2) round(xDatainit{o}(i)-avgshift(f))+win]);
                        if isempty(x)
                            error('ignore point')
                        end
                        xzoomprofile=dataspectra(x);
                        
                        [~,ind]=max(xzoomprofile);
                        ind=ind-win-1;
                        x=max([1 round(xDatainit{o}(i)-avgshift(f))-win+ind]):min([size(spectra,2) round(xDatainit{o}(i)-avgshift(f))+win+ind]);
                        xzoomprofile=dataspectra(x);
                        xzoomprofile=xzoomprofile-xzoomprofile(1);
                        
                        [cf_, b, c, d, a, gof] = fit_gauss(x',xzoomprofile');
                        xfit(i)=b(a==max(a));
                        
                        if gof.rsquare<0.80
                            usemask(i)=false;
                            %                                         figure(1);clf
                            %                                         plot(cf_,x,xzoomprofile,'-x')
                            %                                         hold all
                            %                                         plot(xfit(i),a(a==max(a)),'*')
                            %                                         hold off
                            %                         cf_
                            %                         gof
                            %                         pause;
                        end
                    catch err
                        err
                        usemask(i)=false;
                    end
                end
                if sum(usemask)<4; disp(o); disp(f); error(' '); end
                xDatafit(o,f)={xfit(usemask)};
                xReffibre(o,f)={xRef{o}(usemask)};
                [p(f,:,o),S(f,:,o),mu(f,:,o)] = polyfit(xDatafit{o,f},xReffibre{o,f},3);
            end
        end
        
        save([filename '-autoFittedWave.mat'],'p','S','mu','shift','avgshift')
    end
    return
    %%
    
    flatpayload=load('flat1-001-trace.mat');
    
    flatpayload.flatBlaze(isnan(flatpayload.flatBlaze))=1;
    flatpayload.P2PVariationValues(isnan(flatpayload.P2PVariationValues))=1;
    
    spectra=(fitsread([filename '.fits']))./flatpayload.flatBlaze./flatpayload.P2PVariationValues;
    spectraVar=(fitsread([filename '.fits'],'image',1))./flatpayload.flatBlaze./flatpayload.P2PVariationValues;
    
    spectra(spectra<0)=0;
    %spectraVar(spectraVar<0)=1;
    spectraVar=abs(spectraVar);
    
    
    wavelinear=[];
    for o=1:numOrders;
        for f=1:size(spectra,1);
            wavefit(f,:,o)=polyval(p(f,:,o),1:size(spectra,2),S(f,:,o),mu(f,:,o));
        end
        wavelinear(:,o)=linspace(min(wavefit(10,:,o)),max(wavefit(10,:,o)),size(spectra,2));
    end
    
    
    
    
    speclinear=[];
    for o=1:numOrders;
        for f=1:size(spectra,1);
            speclinear(f,:,o)=interp1(wavefit(f,:,o),spectra(f,:,o),wavelinear(:,o),'pchip',NaN);
            spectraVarlinear(f,:,o)=interp1(wavefit(f,:,o),spectraVar(f,:,o),wavelinear(:,o),'pchip',NaN);
        end
    end
    
    % clip edges
    for o=1:numOrders
        wavelinear_clip{o}=wavelinear(~isnan(sum(speclinear(:,:,o),1)),o);
        speclinear_clip{o}=speclinear(:,~isnan(sum(speclinear(:,:,o),1)),o);
    end
    
    % make single spec for each fibre
    longSpec=[];
    longWave=[];
    for o=1:numOrders
        %longSpec=[longSpec speclinear_clip{o}];
        longWave=[longWave wavelinear_clip{o}'];
    end
    
    longwavelinear=linspace(min(wavefit(:)),max(wavefit(:)),2338*10);
    for o=1:numOrders
        for f=1:size(spectra,1);
            speclinearlong(f,:,o)=interp1(wavefit(f,:,o),spectra(f,:,o),longwavelinear,'pchip',NaN);
            spectraVarlinearlong(f,:,o)=interp1(wavefit(f,:,o),spectraVar(f,:,o),longwavelinear,'pchip',NaN);
        end
    end
    
    finalspeclong=nanmean(speclinearlong,3)';
    finalspecVarlong=nanmean(spectraVarlinearlong,3)';
    
    toclip=isnan(sum(finalspeclong,2));
    
    longwavelinear_clipped=longwavelinear(~toclip);
    finalspecVarlong_clipped=finalspecVarlong(~toclip,:);
    finalspeclong_clipped=finalspeclong(~toclip,:);
    
    header=fitsheader([filename(1:end-11) '.fit']);
    header.IMAGETYP='SPECTRUM';
    header.CRPIX1=round(length(longwavelinear_clipped)/2);
    header.CRVAL1=longwavelinear_clipped(header.CRPIX1);
    header.CTYPE1='Wavelength';
    header.CUNIT1='nm';
    header.CDELT1=mean(diff(longwavelinear_clipped));
    headercell=fitstructure2cell(header);
    
    fitswrite(squeeze(finalspeclong_clipped'),[filename(1:end-11) '-IndivCalSpec.fit'],headercell(8:end,:))
    fitswrite(squeeze(finalspecVarlong_clipped'),[filename(1:end-11) '-IndivCalSpec.fit'],'writemode','append')
    
    fitswrite(squeeze(sum(finalspeclong_clipped,2)'),[filename(1:end-11) '-CombCalSpec.fit'],headercell(8:end,:))
    fitswrite(squeeze(sum(finalspecVarlong_clipped,2)'),[filename(1:end-11) '-CombCalSpec.fit'],headercell(8:end,:))
    
    
    
    
    return
    %plot(wavefit(:,:,1)', spectra(:,:,1)')
    %%
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
    
    
    
    
    
    
    
    
    %%
    subplot(1,2,1)
    plot(longwavelinear_clipped,sum(finalspeclong_clipped,2)./sqrt(sum(finalspecVarlong_clipped,2)))
    subplot(1,2,2)
    plot(longwavelinear_clipped,finalspeclong_clipped./sqrt(finalspecVarlong_clipped))
    
    
    
    
    
    
    
    
    
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
