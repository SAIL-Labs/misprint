%% init reduction class
clc; clear all; set(0,'DefaultFigureWindowStyle','docked')
%set(0,'DefaultFigureWindowStyle','normal')

%for i=1:46
clear s2r
s2r = misprint('AgiDL_Full_band','reference','',... %whyd_1to8
    'plotAlot',true,'numOfOrders',1,'numOfFibers',29,...
    'forceTrace',false,'forceExtract',true,...
    'forceDefineMaskEdge',false,'needsClip',false,...
    'peakcut',0.07,'minPeakSeperation',3,...
    'usecurrentfolderonly',true,...
    'numTraceCol',40,'firstCol',140,'lastCol',300,...
    'parallel',false);%,...
    %'wavesolution','wavefit-ThArg-hgar-preinstall-240-better.mat');%

    
s2r.imdata=s2r.imdata-110;
s2r.imdata(s2r.imdata<0)=NaN;
s2r.filterBadPixels(4,4000,false)
%s2r.removeIntensityGradientInImdata(140)


%%


    
    
%% do stuff to raw image by changing s2r.imdata

%s2r.runDefaultExtraction
s2r.getMaskForIncompleteOrders;
s2r.traceSpectra;

s2r.specCenters=s2r.specCenters-1;

s2r.extractSpectra;
s2r.getP2PVariationsAndBlaze(true);

s2r.plotSpectraFor(1,true,false)
%s2r.plotSingleFibre(1,false,false)


imagesc(s2r.spectraValues)

figure(2)
subplot(1,4,1)
imagesc(s2r.spectraValues)
subplot(1,4,2)
imagesc(s2r.spectraValues./s2r.flatBlaze)
subplot(1,4,3)
imagesc(s2r.spectraValues./s2r.P2PVariationValues)
subplot(1,4,4)
imagesc(s2r.spectraValues./s2r.P2PVariationValues./s2r.flatBlaze)

%%
figure(3)
flatspec=s2r.spectraValues./s2r.P2PVariationValues./s2r.flatBlaze;

for offset=0;%-6:0.1:6
    newflat=[];
    for i=1:29
        newflat(i,:)=interp1([1:320]+(i-1)*offset,flatspec(i,:),-120:1:440);
    end
    
%     subplot(4,1,1:3)
%     imagesc(newflat)
%     
%     limits=[find(isnan(newflat(1,:))==0,1,'first') size(newflat,2)-find(fliplr(isnan(newflat(1,:))==0),1,'first');...
%         find(isnan(newflat(end,:))==0,1,'first') size(newflat,2)-find(fliplr(isnan(newflat(end,:))==0),1,'first')];
%     
%     xlim([min2(limits) max2(limits)])
%     
%     subplot(4,1,4)
%     plot(sum(newflat,1))
%     xlim([min2(limits) max2(limits)])
%     pause
end
%%
newflat(isnan(newflat))=0;
offset=0;
for i=2:29
   [c,lags]=xcorr(newflat(i-1,:),newflat(i,:),20,'coeff');
   [~,ind]=max(c);
   lags(ind)
   
   offset(i)=[offset(i-1) + lags(ind)];
end

    newflat=[];
    for i=1:29
        newflat(i,:)=interp1([1:320]+(i-1)*offset(i),flatspec(i,:),-120:1:440);
    end
    
    subplot(4,1,1:3)
    imagesc(newflat)
    
    limits=[find(isnan(newflat(1,:))==0,1,'first') size(newflat,2)-find(fliplr(isnan(newflat(1,:))==0),1,'first');...
        find(isnan(newflat(end,:))==0,1,'first') size(newflat,2)-find(fliplr(isnan(newflat(end,:))==0),1,'first')];
    
    %xlim([min2(limits) max2(limits)])
    
    subplot(4,1,4)
    plot(sum(newflat,1))
    %xlim([min2(limits) max2(limits)])


