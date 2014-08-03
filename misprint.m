classdef misprint < handleAllHidden
    %MISPRINT MultIorder SPectroscopic ReductIoN Tool
    %
    % Properties are set using key value pairs.
    %
    % Sample Usage:
    %     s2r = misprint('sciencespectrum','reference','flatspectrum','plotAlot',true,...
    %         'usecurrentfolderonly',true,...
    %         'numOfOrders',14,'numOfFibers',29,...
    %         'forceTrace',false,'forceExtract',false,...
    %         'forceDefineMaskEdge',false,'needsMask',false,...
    %         'peakcut',0.07,'minPeakSeperation',3,...
    %         'numTraceCol',40,'firstCol',140,'lastCol',300,...
    %         'parallel',false);
    %
    %     self.getMaskForIncompleteOrders;
    %     self.traceSpectra;
    %     self.extractSpectra;
    %     self.getP2PVariationsAndBlaze(false);
    %
    %     self.plotSpectraFor(1:14,true,false)
    %
    % Copyright (C) Chris Betters 2012-2014
    
    properties
        targetBaseFilename, % base filename of target (file with spectra to be extracted) fits file.
        targetPath, % path to target fits file.
        rootDirectory, % path to current directory, should equal [pwd '/'].
        referenceBaseFilename, % base filename of reference fits file (i.e. a flat frame).
        referencePath, % path to reference fits file.
        
        targetHeader, % structure with target fits header.
        referenceHeader, % structure with reference fits header.
        
        SpectraFitsSaveFileName, % filename of fits file to save extracted spectra to.
        ReferenceSpectraFitsSaveFileName, % spectra previously extracted from the reference fits file.
        FlatReferenceSpectraFitsSaveFileName, % filename of fits file to for use as a flat reference, should be just pixel to pixel variations.
        
        spectraTracePath, % path to previously saved trace data for reference fits file.
        
        useReference, % flag to indicate if valid reference data has been set/ffound.
        plotAlot, % flag to show raw image and plots during tracing. It can plot alot.
        forceTrace, % flag to force a trace of spectra in the current image. Can not be set with a reference.
        forceExtract, % flag to force a new extraction of the current image. This will overwrite the default save file if it exists.
        forceDefineMaskEdge, % flag to force a new definition of the mask/clipping region.
        needsMask, % flag to indicate if image requires clipping.
        clipping, % vector of pixels from [left top right bottom] to clip.
        parallel, % flag tin indicate if parallel compuyting toolbox should/can be used.
        minPeakSeperation, % min peak seperation for tracing and peakfinder
        
        numOfOrders, % number of diffraction orders in image.
        numOfFibers, % number of spectra (fibres) in each order.
        
        gain, % gain (e-/adu) read from fits file
        readNoise, % read noise (rms e-) read from fits file
        dispAxis, % axis of primary dispersion read from fits file
        
        imdata, % target image data.
        imvariance, % estiamted variance for target image data.
        mask, % mask of clipped regions.
        imdim, % size of imdata (equals size(imdata))
        
        spectraValues, % extracted spectra values.
        spectraVar, % var for extracted spec values.
        backgroundValues, % background value from extraction
        
        referenceSpectraValues, % extracted spectraValues of reference/flat
        P2PVariationValues, % pixel to pixel variation from reference.
        flatBlaze, % estimated blaze from reference.
        
        wavefit,% matfile with wavefit
        diffractionOrder, % estimated diffraction order from wavefit
        
        numTraceCol, % number of columns to use when tracing.
        firstCol, % first column of trace
        lastCol, % last column of trace
        
        orderEdges,  % detected edges of the orders
        specCenters, % polynoimail interpolated fitted y axis centeres of the spectra from gaussain fit.
        specWidth,   % polynomial interpolated width of the spectra from gaussain fit.
        meanSpecWidth, % mean width of spceetra in each order.
        meanOrderWidth, % mean width of each order
        
        fittedCenters, % center of spectrum (vertical) from gaussian fit
        fittedCol, % column used to get progile for fit
        fittedWidth, % width of specturm (vertical) from gaussain fit
        fittedParamters, % all fit paramaters from trace
        
        usecurrentfolderonly, % flag to note use my maxiumDL/PIMMS echelle file structure.
        peakcut, % MINPEAKHEIGHT for spectra tracing detection. (fraction of mean of current profile).
    end
    
    methods
        function self=misprint(targetBaseFilename,varargin)
            % init the MISPRINT class. Pasre all inputs, load relevant files.
            
            %% parse inputs
            p = inputParser;
            
            p.addRequired('targetBaseFilename', @(x) ischar(x));
            p.addParamValue('reference','');
            p.addParamValue('forceTrace'            , false, @(x) islogical(x));
            p.addParamValue('forceExtract'          , false, @(x) islogical(x));
            p.addParamValue('plotAlot'      , false, @(x) islogical(x));
            p.addParamValue('forceDefineMaskEdge', false, @(x) islogical(x));
            p.addParamValue('needsMask', true, @(x) islogical(x));
            p.addParamValue('numOfOrders', 15, @(x) isnumeric(x));
            p.addParamValue('numOfFibers', 19, @(x) isnumeric(x));
            p.addParamValue('usecurrentfolderonly',false, @(x) islogical(x));
            p.addParamValue('peakcut', 0.8, @(x) isnumeric(x));
            p.addParamValue('parallel',true, @(x) islogical(x));
            p.addParamValue('numTraceCol', 10, @(x) isnumeric(x));
            p.addParamValue('dispAxis', [], @(x) isnumeric(x));
            p.addParamValue('wavesolution', '', @(x) ischar(x));
            p.addParamValue('minPeakSeperation', 3, @(x) isnumeric(x));
            p.addParamValue('firstCol', 0, @(x) isnumeric(x));
            p.addParamValue('lastCol', 0, @(x) isnumeric(x));
            p.addParamValue('clipping',[0 0 0 0], @(x) isnumeric(x) && length(x)==4)
            
            
            p.parse(targetBaseFilename,varargin{:});
            
            self.numOfOrders=p.Results.numOfOrders;
            self.numOfFibers=p.Results.numOfFibers;
            self.forceExtract=p.Results.forceExtract;
            self.forceTrace=p.Results.forceTrace;
            self.plotAlot=p.Results.plotAlot;
            self.forceDefineMaskEdge=p.Results.forceDefineMaskEdge;
            self.needsMask=p.Results.needsMask;
            
            self.usecurrentfolderonly=p.Results.usecurrentfolderonly;
            self.peakcut=p.Results.peakcut;
            
            self.parallel=p.Results.parallel;
            self.numTraceCol=p.Results.numTraceCol;
            self.firstCol=p.Results.firstCol; %if zero set to 20 minus miage size (at end)
            self.lastCol=p.Results.lastCol; %if zero set to 20 minus miage size (at end of init)
            
            self.dispAxis=p.Results.dispAxis;
            
            self.minPeakSeperation=p.Results.minPeakSeperation;
            
            self.clipping=p.Results.clipping;
            
            if ~isempty(p.Results.wavesolution)
                matpayload=load(p.Results.wavesolution);
                self.wavefit=matpayload.wavefit;
                self.diffractionOrder=round(2*1e-3/31.6*cosd(5)*sind(63.2)./(mean(squeeze(mean(self.wavefit,2)),1)*1e-10)');
            end
            
            %% start matlabpool if parallel computing tool box avaliable
            if self.parallel
                if license('test', 'distrib_computing_toolbox')
                    if ~matlabpool('size')
                        matlabpool(4)
                    end
                else
                    warning('MISPRINT:init:useDistribComputingToolbox:notAvalaible','The parallel computing toolbox is not avaliable, but has been requested.')
                end
            end
            
            %% root path
            self.rootDirectory=[pwd '/'];
            
            %% main reduction target path construction
            self.targetBaseFilename=p.Results.targetBaseFilename;
            
            if ~self.usecurrentfolderonly
                self.targetPath = [self.rootDirectory self.targetBaseFilename '/reduced/' self.targetBaseFilename '-master.fit'];
            else
                self.targetPath = [self.rootDirectory self.targetBaseFilename '.fit'];
            end
            % check the file is a valid fits.
            self.checkForReducedFitsAt(self.targetPath);
            
            % get fits header
            self.targetHeader=fitsheader(self.targetPath);
            
            %% reference target path construction
            self.referenceBaseFilename=p.Results.reference;
            
            if isempty(self.referenceBaseFilename)
                self.useReference     = false;
                self.spectraTracePath = [self.rootDirectory self.targetBaseFilename '-trace.mat'];
            else
                self.useReference     = true;
                self.spectraTracePath = [self.rootDirectory self.referenceBaseFilename '-trace.mat'];
                if ~self.usecurrentfolderonly
                    self.referencePath    = [self.rootDirectory self.referenceBaseFilename '/reduced/' self.referenceBaseFilename '-master.fit'];
                else
                    self.referencePath    = [self.rootDirectory self.referenceBaseFilename '.fit'];
                end
                % check the file is a valid fits.
                self.checkForReducedFitsAt(self.targetPath);
                
                % get fits header
                self.referenceHeader  = fitsheader(self.referencePath);
                
                
                self.ReferenceSpectraFitsSaveFileName=[self.referenceBaseFilename '-1D-spectra.fits'];
                self.referenceSpectraValues=fitsread(self.ReferenceSpectraFitsSaveFileName);
                
                self.FlatReferenceSpectraFitsSaveFileName=[self.referenceBaseFilename '-flattend-1D-spectra.fits'];
                
                assertWarn(isfield(self.referenceHeader,'IMAGETYP') && strcmp(self.referenceHeader.IMAGETYP,'Flat Frame'),...
                    'MISPRINT:init:referenceNotAFlat',...
                    'Reference Frame has not been tagged as a flat in fits header')
            end
            
            %% 1D spectra filenames
            self.SpectraFitsSaveFileName=[self.targetBaseFilename '-1D-spectra.fits'];
            
            %% check for required cards in fits header and read the values. add defaults where unavaliable.
            if isfield(self.targetHeader,'READNOIS')
                self.readNoise=self.targetHeader.READNOIS;
            elseif isfield(self.targetHeader,'RO_NOISE')
                self.readNoise=self.targetHeader.RO_NOISE;
            else
                self.readNoise=10; % atik default
                %                 if strcmp(self.targetHeader.INSTRUME,'ArtemisHSC')
                %                     fitsAddHeaderKeyword(self.targetPath,'READNOIS',self.readNoise,' ')
                %                 end
            end
            
            if isfield(self.targetHeader,'GAIN')
                self.gain=self.targetHeader.GAIN;
            elseif isfield(self.targetHeader,'RO_GAIN')
                self.gain=self.targetHeader.RO_GAIN;
            else
                self.gain=0.46; % atik default
                %                 if strcmp(self.targetHeader.INSTRUME,'ArtemisHSC')
                %                     fitsAddHeaderKeyword(self.targetPath,'GAIN',self.gain,' ')
                %                 end
            end
            
            if isempty(self.dispAxis)
                if isfield(self.targetHeader,'DISPAXIS')
                    self.dispAxis=self.targetHeader.DISPAXIS;
                else
                    self.dispAxis=1; % atik default
                    %if strcmp(self.targetHeader.INSTRUME,'ArtemisHSC')
                    %    fitsAddHeaderKeyword(self.targetPath,'DISPAXIS',self.dispAxis,' ');
                    %end
                end
            end
            %% load misprint, and orinateate so echelle dispersion is horizontal
            self.imdata=fitsread(self.targetPath);
            self.dispAxis
            if self.dispAxis==1
                self.imdata=(self.imdata'); %fliplr
            end
            
            if sum(self.clipping)
                %[left top right bottom]
                self.imdata=self.imdata(max([1 self.clipping(2)]):end-self.clipping(4),max([1 self.clipping(1)]):end-self.clipping(3));
                if ~isempty(self.wavefit)
                    self.wavefit=self.wavefit(:,max([1 self.clipping(1)]):end-self.clipping(3),:);
                end
            end
            
            %self.imdata=rot90(self.imdata,2);
            self.imdim=size(self.imdata);
            self.imvariance=(self.readNoise/self.gain)^2 + abs(self.imdata) / self.gain; %http://cxc.cfa.harvard.edu/mta/ASPECT/aca_read_noise/
            
            %% trace col
            if ~self.lastCol
                self.lastCol=self.imdim(2)-20;
            end
            if ~self.firstCol
                self.firstCol=20;
            end
            
        end
        
        function runDefaultExtraction(self)
            % run default set of extraction commands
            
            self.getMaskForIncompleteOrders;
            self.traceSpectra;
            self.extractSpectra;
            self.getP2PVariationsAndBlaze
        end
        
        function spectraValues=extractSpectra(self)
            % extract spectra using trace - each order done individually (faster). 
            
            if ~exist(self.SpectraFitsSaveFileName,'file') || self.forceExtract
                spectraValues=zeros(self.numOfFibers,self.imdim(2),self.numOfOrders);
                spectraVar=zeros(self.numOfFibers,self.imdim(2),self.numOfOrders);
                warning('spec width is fixed')
                assertWarn(self.forceExtract,...
                    'MISPRINT:extractSpectra:forceExtractFlagSet',...
                    'Force extraction flag set, starting extraction. Data will be overwritten')
                
                for order=1:self.numOfOrders
                    spectra=zeros(self.numOfFibers,self.imdim(2));
                    specVar=zeros(self.numOfFibers,self.imdim(2));
                    
                    disp(['Extracting Order: ' num2str(order)])
                    for col=1:self.imdim(2)
                        orderCenter=mean(self.specCenters(order,:,col));
                        
                        profileApeture=max(round(orderCenter-self.meanOrderWidth/2),1):...
                            min(round(orderCenter+self.meanOrderWidth/2),self.imdim(1));
                        
                        orderProfile=self.imdata(profileApeture,col);
                        varProfile=self.imvariance(profileApeture,col)';
                        
                        %         x=repmat(profileApeture',[19,1]);
                        %         bsxfun(@minus,x,squeeze(specCenters(order,:,col)))
                        [spectra(:,col), specVar(:,col)]=self.MPDoptimalExt(profileApeture,orderProfile,varProfile,...
                            squeeze(self.specCenters(order,:,col)),squeeze(self.specWidth(order,:,col)),...
                            self.readNoise/self.gain);
                    end
                    spectraValues(:,:,order)=spectra;
                    spectraVar(:,:,order)=specVar;
                end
                
                spectra1DHDR=createcards('NUMORDER',self.numOfOrders,'number of orders');
                spectra1DHDR.addcard('NUMFIBER',self.numOfFibers,'number of fibers')
                spectra1DHDR.addcard('TRACE',self.spectraTracePath,' ')
                
                fitswrite(spectraValues,self.SpectraFitsSaveFileName,spectra1DHDR.cards)
                fitswrite(spectraVar,self.SpectraFitsSaveFileName,'writemode','append')
            else
                disp(['Pre-existing extraction data found at: ' self.SpectraFitsSaveFileName])
                
                spectraValues=fitsread(self.SpectraFitsSaveFileName);
                spectraVar=fitsread(self.SpectraFitsSaveFileName,'image');
            end
            self.spectraValues=spectraValues;
            self.spectraVar=spectraVar;
        end
        
        function traceSpectra(self,varargin)
            % trace spectra from flat.
            %
            % optional inputs misprint.traceSpectra(inputimage,numOfOrders,numOfFibers)
            % inputimage is same format as misprint.imdata
            
            %% inital setup
            % if preexisting trace exists it is loaded (unless forceTrace set)
            if nargin==1
                if (~exist(self.spectraTracePath,'file') || self.forceTrace ) && ~self.useReference
                    
                    assertWarn(self.forceTrace & exist(self.spectraTracePath,'file'),...
                        'MISPRINT:traceSpectra:TraceForced',...
                        'Tracing was forced, this will overwrite previous trace.')
                    
                    if ~exist(self.spectraTracePath,'file'); disp(['Tracefile: ' self.spectraTracePath ' does not exist.']);end
                else
                    assert(~(self.forceTrace & self.useReference),...
                        'MISPRINT:traceSpectra:TraceForcedWithUseReferenceSet',...
                        'Tracing can not be forced when useReference is set')
                    
                    assert(~(self.useReference & ~exist(self.spectraTracePath,'file')),...
                        'MISPRINT:traceSpectra:ReferecnceTraceFileNotFound',...
                        [self.spectraTracePath ' was not found and is required as a reference. Aborting extraction.'])
                    
                    load(self.spectraTracePath,'specCenters','specWidth','orderWidth','orderEdges','means','columns','widths','fitxs')
                    
                    if self.useReference
                        disp(['Using reference trace: ' self.spectraTracePath])
                    else
                        disp(['Using previous trace: ' self.spectraTracePath])
                    end
                    
                    self.meanSpecWidth=squeeze(mean(specWidth,2));
                    self.meanOrderWidth=squeeze(mean(orderWidth,2));
                    self.specCenters=specCenters;
                    self.specWidth=specWidth;
                    self.orderEdges=orderEdges;
                    
                    self.fittedCenters=means;
                    self.fittedCol=columns;
                    self.fittedWidth=widths;
                    self.fittedParamters=fitxs;
                    
                    return % end function call after loading data
                end
            end
            
            %% load data into local variables
            x=1:self.imdim(1);
            if nargin==4
                inputimage = varargin{1};
                numOfOrders = varargin{2};
                numOfFibers = varargin{3};
            else
                inputimage=self.imdata;
                numOfOrders = self.numOfOrders;
                numOfFibers = self.numOfFibers;
            end
            
            imdata=inputimage.*self.mask;
            
            
            %% find orders
            if (self.numTraceCol>=self.imdim(2))
                columns=1:self.imdim(2);
                warning('MISPRINT:fitAllOfTheThings','You just asked for a fit to every column.')
                reply = input('Are your sure?? Y/N [Y]: ', 's');
                if ~strcmpi('Y',reply)
                    error('MISPRINT:traceSpectra:userInterpupt','MISPRINT termintated in traceSpectra.')
                end
            else
                columns=round(linspace(self.firstCol, self.lastCol, self.numTraceCol));
            end
            
            imcol=imdata(:,columns); % sliced image
            
            disp('Running order tracer. This may take some time.')
            for i=1:length(columns)
                [yp,index]=findpeaks(imcol(:,i),'NPEAKS',numOfOrders*numOfFibers,'MINPEAKHEIGHT',max(imcol(:,i))*self.peakcut,'MINPEAKDISTANCE',self.minPeakSeperation);
                if self.plotAlot
                    figure(i);clf
                    plot(x,imcol(:,i),index,yp,'xr');
                    line([1 length(imcol(:,i))],[max(imcol(:,i)) max(imcol(:,i))]*self.peakcut)
                    title([num2str(columns(i))])
                end
                if numOfOrders==1
                    orderWidth=self.imdim(1);
                    orderCenter=round(self.imdim(1)/2);
                    orderEdges(:,i)=[1 self.imdim(1)];
                else
                    orderWidth=diff(index(1:numOfFibers:end));
                    orderCenter=mean([index(numOfFibers:numOfFibers:end) index(1:numOfFibers:end)],2);
                    %error(' ')
                    orderEdges(:,i)=[orderCenter(1)-orderWidth(1)/2;...
                        mean([index(numOfFibers:numOfFibers:end-numOfFibers)...
                        index(numOfFibers+1:numOfFibers:end)],2); orderCenter(end)+orderWidth(end)/2];
                end
                
                %                 if self.plotAlot
                %                     plot(imcol(:,i))
                %                     hold on
                %                     %line(repmat(orderEdges(:,i)',[2,1]),[zeros(1,size(orderEdges(:,i),2)); ones(1,size(orderEdges(:,i),2))*max(imcol)])
                %                     hold off
                %                 end
            end
            
            %% trace orders
            % fit gaussian to profile in columns for each order.
            fitxs=zeros(numOfOrders,3*numOfFibers+1,length(columns));
            for i=1:length(columns)
                for order=1:numOfOrders
                    disp(['Column:' num2str(columns(i)) ' | Fitting Spectra in Order: ' num2str(order)])
                    
                    orderProfileX=round(max(orderEdges(order,i),1):min(orderEdges(order+1,i),self.imdim(1)))';
                    orderProfile=imcol(orderProfileX,i);
                    orderProfile=orderProfile/max(orderProfile);
                    
                    
                    %                     [~, means(order,:,i), widths(order,:,i), fitxs(order,:,i)] = ...
                    %                         fitNGaussainsAlt(numOfFibers,orderProfileX, orderProfile,self.peakcut);
                    
                    [~, means(order,:,i), widths(order,:,i), fitxs(order,:,i)] = ...
                        fitNGaussains(numOfFibers,orderProfileX, orderProfile,self.peakcut);
                    
                    if self.plotAlot
                        figure(i);clf;
                        %subplot(5,4,columns)
                        plot(orderProfileX,sum(nGausFunc(fitxs(order,:,i),orderProfileX,numOfFibers),2),'r-',...
                            orderProfileX,orderProfile,'-')
                        title(['Order: ' num2str(order) ' Column: ' num2str(columns(i))])
                        %pause(0.1)
                    end
                end
            end
            
            specCenters=polyfitwork(self.imdim,means,columns,2,0);
            specWidth=polyfitwork(self.imdim,widths,columns,2,0);
            meanSpecWidth=squeeze(mean(self.specWidth,2));
            
            self.fittedCenters=means;
            self.fittedCol=columns;
            self.fittedWidth=widths;
            
            self.meanSpecWidth=meanSpecWidth;
            self.specCenters=specCenters;
            self.specWidth=specWidth;
            self.orderEdges=orderEdges;
            self.fittedParamters=fitxs;
            self.meanOrderWidth=squeeze(mean(orderWidth,2));
            
            save(self.spectraTracePath,'specCenters','specWidth','orderWidth','orderEdges','means','columns','widths','fitxs')
        end
        
        function getMaskForIncompleteOrders(self)
            %  get mask for incomplete orders
            
            if ~self.needsMask
                self.mask=ones(self.imdim);
                return % no clipe, so mask is ones.
            end
            
            if isfield(self.targetHeader,'CLIPTL') && isfield(self.targetHeader,'CLIPTR') && isfield(self.targetHeader,'CLIPBL') && isfield(self.targetHeader,'CLIPBR')
                topEdges=[self.targetHeader.CLIPTL self.targetHeader.CLIPTR];
                bottomEdges=[self.targetHeader.CLIPBL self.targetHeader.CLIPBR];
            else
                self.forceDefineMaskEdge=true; % overide default as clipping is needed, and data not defined
            end
            
            if self.forceDefineMaskEdge
                echfig=figure(1);
                imagesc(self.imdata);
                axis([1 self.imdim(2) 1 self.imdim(1)*0.3]) % show top third of image
                [~,y]=getpts(echfig);
                topEdges=[y(1) y(2)];
                axis([1 self.imdim(2) self.imdim(1)-self.imdim(1)*0.3 self.imdim(1)]) % show bottom third of image
                [~,y]=getpts(echfig);
                bottomEdges=[y(1) y(2)];
            end
            
            % make mask of image to exclude incomplete orders
            xi=[0; self.imdim(2);              self.imdim(2);                  0];
            yi=[0;        0;    topEdges(2); topEdges(1)];
            BW1 = roipoly(self.imdata,xi,yi);
            
            xi=[       0; self.imdim(2);           self.imdim(2);                  0];
            yi=[self.imdim(1); self.imdim(1);     bottomEdges(2);     bottomEdges(1)];
            
            BW2 = roipoly(self.imdata,xi,yi);
            self.mask=~BW1 & ~BW2;
            
            if self.forceDefineMaskEdge
                imagesc(self.imdata.*self.mask)
                reply = input('Should I add to Fits Header Y/N [N]: ', 's');
                if isempty(reply)
                    reply = 'N';
                end
                if strncmpi(reply,'Y',1)
                    disp('saving clips to header')
                    fitsAddHeaderKeyword(self.targetPath,'CLIPTL',topEdges(1),' ');
                    fitsAddHeaderKeyword(self.targetPath,'CLIPTR',topEdges(2),' ');
                    fitsAddHeaderKeyword(self.targetPath,'CLIPBL',bottomEdges(1),' ');
                    fitsAddHeaderKeyword(self.targetPath,'CLIPBR',bottomEdges(2),' ');
                end
            end
            
            if self.plotAlot
                figure(1)
                imagesc(self.imdata.*self.mask)
            end
            
        end
        
        function getP2PVariationsAndBlaze(self,varargin)
            % get smoothed version of flat spectrum (ie blaze) and pixel to
            % pixel variations (flatspectrum./smooth flat spectrum)
            %
            %  load reference
            if length(varargin)==2
                referenceFile=varargin{2};
            elseif self.useReference
                referenceFile=self.ReferenceSpectraFitsSaveFileName;
            end
            
            if ~isempty(varargin)
                force=varargin{1};
            else
                force=false;
            end
            
            matpayload=load(self.spectraTracePath,'flatBlaze','P2PVariationValues');
            if ~force && isfield(matpayload,'flatBlaze') && isfield(matpayload,'P2PVariationValues')
                self.flatBlaze=matpayload.flatBlaze;
                self.P2PVariationValues=matpayload.P2PVariationValues;
            else
                assertWarn(force,'MISPRINT:getP2PVariationsAndBlaze:forced','P2P and blaze forced')
                mask=ones(self.numOfFibers,self.imdim(2));
                mask(end-50:end)=NaN;
                mask(1:50)=NaN;
                
                spectraValues=self.spectraValues;
                for or=1:self.numOfOrders
                    spectraValues(:,:,or)=bsxfun(@rdivide,spectraValues(:,:,or).*mask,mean(spectraValues(:,:,or)')');
                end
                flatBlaze=zeros(size(self.spectraValues));
                P2PVariationValues=zeros(size(self.spectraValues));
                
                for f=1:self.numOfFibers
                    for or=1:self.numOfOrders
                        flatBlaze(f,:,or)=csaps(1:self.imdim(2),spectraValues(f,:,or),1e-9,1:self.imdim(2));
                        %flatBlaze(f,:,or)=smooth(spectraValues(f,:,or),100);
                    end
                end
                P2PVariationValues=spectraValues./flatBlaze;
                self.flatBlaze=flatBlaze;
                self.P2PVariationValues=P2PVariationValues;
                save(self.spectraTracePath,'flatBlaze','P2PVariationValues','-append')
                %error(' ')
            end
        end
        
        function plotSpectraFor(self,orders,shouldFlat,shouldP2PV)
            % plot spectra orders specifed. three arguments orders,shouldFlat,shouldP2PV
            
            if shouldFlat && shouldP2PV
                FlattenedSpectra=self.spectraValues./self.flatBlaze./self.P2PVariationValues;
            elseif shouldFlat && ~shouldP2PV
                FlattenedSpectra=self.spectraValues./self.flatBlaze;
            elseif ~shouldFlat && shouldP2PV
                FlattenedSpectra=self.spectraValues./self.P2PVariationValues;
            else
                FlattenedSpectra=self.spectraValues;
            end
            
            if isempty(self.wavefit)
                for order=orders
                    figure(order);clf;
                    subplot(1,2,1)
                    %imagesc(log10(self.imdata-min2(self.imdata)+1))
                    imagesc(self.imdata)
                    %set(gca, 'CLim',[0 1000])
                    ylim([min2(squeeze(self.specCenters(order,:,:)))-50 max2(squeeze(self.specCenters(order,:,:)))+50])
                    hold on
                    plot(1:self.imdim(2),squeeze(self.specCenters(order,:,:)),'k','LineWidth',1.5)
                    %axis image
                    
                    subplot(1,2,2)
                    for f=1:self.numOfFibers
                        FlattenedSpectraNorm(f,:,order)=FlattenedSpectra(f,:,order)/(max(FlattenedSpectra(f,:,order)));
                        plot(1:self.imdim(2),FlattenedSpectraNorm(f,:,order)+(self.numOfFibers-f)*1)
                        hold all
                    end
                end
                hold off
            else
                
                for order=orders
                    figure(order);clf;
                    
                    subplot(1,2,1)
                    imagesc(log10(self.imdata-min2(self.imdata)+1))
                    %imagesc(self.imdata)
                    %set(gca, 'CLim',[0 1000])
                    ylim([min2(squeeze(self.specCenters(order,:,:)))-50 max2(squeeze(self.specCenters(order,:,:)))+50])
                    hold on
                    plot(1:self.imdim(2),squeeze(self.specCenters(order,:,:)),'k','LineWidth',1.5)
                    title(['Order ' num2str(self.diffractionOrder(order))])
                    xlabel('Primary-dispersion axis (pixels)')
                    ylabel('Cross-dispersion axis (pixels)')
                    
                    subplot(1,2,2)
                    %                     for f=1:self.numOfFibers
                    %                         FlattenedSpectraNorm(f,:,order)=FlattenedSpectra(f,:,order)/max(FlattenedSpectra(f,:,order));
                    %                         plot(self.wavefit(f,:,order),FlattenedSpectraNorm(f,:,order)+(self.numOfFibers-f)*0.25)
                    %                         hold all
                    %                     end
                    
                    plot(self.wavefit(:,:,order)',FlattenedSpectra(:,:,order)')
                    xlabel('Wavelength (nm)')
                end
                hold off
                
            end
        end
        
        function plotSingleFibre(self,f,shouldFlat,shouldP2PV)
            % plot spectra for signle fibre across multiple orders, three arguments fibre,shouldFlat,shouldP2PV
            if shouldFlat && shouldP2PV
                FlattenedSpectra=self.spectraValues./self.flatBlaze./self.P2PVariationValues;
            elseif shouldFlat && ~shouldP2PV
                FlattenedSpectra=self.spectraValues./self.flatBlaze;
            elseif ~shouldFlat && shouldP2PV
                FlattenedSpectra=self.spectraValues./self.P2PVariationValues;
            else
                FlattenedSpectra=self.spectraValues;
            end
            
            if isempty(self.wavefit)
                plot([1:self.imdim(2)],squeeze(FlattenedSpectra(f,:,:)))
            else
                plot(squeeze(self.wavefit(f,:,:)),squeeze(FlattenedSpectra(f,:,:)))
            end
        end
        
        function filterBadPixels(self,Nsigma,thresh,shouldPlot)
            % filter bad pixels. three arguments Nsigma,thresh,shouldPlot
            im=self.imdata;
            im(im<=0)=1;

            imdiff=medfilt2(im,[2 2])./im; % try and highlight odd pixels
            
            imdiff=imdiff-mean2(imdiff); % set mean to zero
            
            bad1=imdiff>thresh; % very larger value can bias std, so clip them.
            
            badpixel=(imdiff>std2(imdiff(~bad1))*Nsigma | imdiff<-std2(imdiff(~bad1))*Nsigma ) | bad1;
            
            self.imdata(badpixel)=NaN;
            
            self.imdata=inpaint_nans(self.imdata);
            
            if shouldPlot
                figure(shouldPlot);clf
                [badx, bady]=find(badpixel);
                
                imagesc(self.imdata)
                hold on
                shouldPlot(bady,badx,'wx')
                hold off
            end
            
        end
        
        function removeIntensityGradientInImdata(self,avgWin)
            % smooth whole image, then divided original by that. Usefull for to
            % improve flat tracing. one arguments avgWin (window for smoothing)
            PSF = fspecial('average', [1 1]*avgWin);
            
            blurred = imfilter(self.imdata, PSF, 'conv', 'symmetric');
            blurred=blurred/mean2(blurred);
            
            if self.plotAlot
                subplot(1,3,1)
                imagesc(self.imdata)
                subplot(1,3,2)
                imagesc(blurred)
                subplot(1,3,3)
                imagesc(self.imdata./blurred)
            end
            self.imdata=self.imdata./blurred;
        end
        
        function getBackgroundBetweenOrders(self)
            self.imdata(self.imdata<0)=0;
            locs=self.orderEdges';
            locs(locs>3362)=3362;
            
            imagesc(log10(self.imdata))
            hold on;
            plot(self.fittedCol,locs','bx')
            hold off;
            pks=[];
            
            
            filtedimdata=medfilt2(self.imdata);
            
            %locs(89,:)=(locs(88,:) + locs(90,:)) / 2;
            
            for i = 1:size(locs,2)
                %    pks(i,:)=self.imdata(self.fittedCol,round(self.orderEdges(,:)));
                p = impixel(filtedimdata,self.fittedCol,locs(:,i)');
                pks(:,i) = p(:,1);
            end
            
            
            
            cols=self.fittedCol;
            
            
            
            %% scattered light estimate
            
            invertedimdata=(1./(self.imdata)).*self.mask;
            invertedimdata(isinf(invertedimdata))=0;
            figure(1)
            imagesc(log10(invertedimdata))
            x=1:self.imdim(1);
            
            
            inverpks=1./pks;
            cols2=repmat(cols,[size(locs,2),1])';
            figure(3);clf
            h(2)=surface(cols2,locs,pks,'EdgeColor','none');
            xlim([1 self.imdim(2)])
            ylim([1 self.imdim(1)])
            set(get(h(2),'Parent'),'YDir','reverse')
            
            
            figure(2);clf
            sfun=scateringTestFit(cols2, locs, inverpks);
            
            
            
            figure(4);clf
            [XI,YI]=meshgrid(1:self.imdim(2), 1:self.imdim(1));
            
            subplot(1,2,2)
            imagesc(1./feval(sfun,XI,YI).*self.mask)
            title('Estimated Scattering (from Inter-Order Regions)')
            
            subplot(1,2,1)
            h(2)=surface(cols2,locs,1./pks,'EdgeColor','none');
            xlim([1 self.imdim(2)])
            ylim([1 self.imdim(1)])
            set(get(h(2),'Parent'),'YDir','reverse')
            
            
            %self.imdata=self.imdata-feval(sfun,X,Y)
            imagesc(self.imdata-1./feval(sfun,XI,YI))
            hold on; plot(cols2,locs,'wx');hold off
            title('PIMMS Echelle Detector Image')
            
            %self.imdata=self.imdata-1./feval(sfun,XI,YI)
            
            
            %%
            
            imagesc(log10(self.imdata))
            
            return
            self.forceTrace=true;
            self.forceExtract=true;
            
            self.getMaskForIncompleteOrders;
            self.traceSpectra;
            %self.specCenters=self.specCenters;
            self.extractSpectra;
            self.getP2PVariationsAndBlaze
            set(0,'DefaultFigureWindowStyle','docked')
            
        end
        
        function [spectraValues, spectraErrors]=MPDoptimalExt(~,dataRows,orderProfile,varProfile,specCenters,specWidth,RN)
            % Multi-Profile Deconvolution Optimal Extraction as described by Sharp & Birchall (2010)
            % 
            % paper: Sharp R., Birchall M. N. (2010) Optimal Extraction of Fibre Optic Spectroscopy. PASA 27, pp. 91-103.
            %        http://dx.doi.org/10.1071/AS08001
            
            if iscolumn(orderProfile); orderProfile=orderProfile'; end
            if iscolumn(dataRows); dataRows=dataRows'; end
            if isrow(specCenters);specCenters=specCenters'; end
            if isrow(specWidth);specWidth=specWidth'; end
            
            if isempty(varProfile)
                varProfile=ones(size(orderProfile));
            end
            
            orderProfile=orderProfile-min(orderProfile);
            
            phi=getPhi(dataRows,specCenters,specWidth,ones(length(specCenters),1));
            
            sigmaweightedPhi=bsxfun(@rdivide,phi,sqrt(varProfile))';
            
            c=phi*sigmaweightedPhi;
            b=(orderProfile*sigmaweightedPhi)';
            
            ce=phi*phi';
            be=((varProfile-RN^2)*phi')';
            
            spectraValues=c\b;
            spectraErrors=ce\be;
            
            function phi=getPhi(dataRows,specCenters,specWidth,specPeaks)
                phi1=bsxfun(@minus,repmat(dataRows,[length(specCenters),1]),specCenters);
                phi2=bsxfun(@rdivide, phi1, specWidth);
                phi3=bsxfun(@times, exp(-(phi2).^2), 1./(specWidth*sqrt(pi)));
                phi=bsxfun(@times, phi3, specPeaks);
                
                phi(phi<1e-8)=0;
            end
        end
    end
    
    methods (Static, Access = private)
        function answer=checkForReducedFitsAt(path)
            % check for a reduced target
            try
                import matlab.io.*
                fptr = fits.openFile(path);
                fits.closeFile(fptr);
                answer=1;
            catch err
                if strcmp(err.identifier,'MATLAB:imagesci:fits:libraryError')
                    error('MISPRINT:checkForReducedTarget:fitsOpenError','Reduced fits file does not exist.')
                else
                    rethrow(err)
                end
            end
        end
    end
end