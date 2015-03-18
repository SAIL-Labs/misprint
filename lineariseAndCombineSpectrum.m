function lineariseAndCombineSpectrum(targetfile,wavematfile,flatmatfile)

%tagetfile - of format <basefilename>-1D-spectra.fits



flatpayload=load(flatmatfile);
flatpayload.flatBlaze(isnan(flatpayload.flatBlaze))=1;
flatpayload.P2PVariationValues(isnan(flatpayload.P2PVariationValues))=1;

spectra=fliplr(fitsread(targetfile))./flatpayload.flatBlaze./flatpayload.P2PVariationValues;
spectraVar=fliplr(fitsread(targetfile,'image',1))./flatpayload.flatBlaze./flatpayload.P2PVariationValues;

for or=1:10
    spectraVar(:,:,or)=bsxfun(@rdivide,spectra(:,:,or),max(spectra(:,:,or)')');
    spectra(:,:,or)   =bsxfun(@rdivide,spectra(:,:,or),max(spectra(:,:,or)')');
end
error(' ')
%wavelinear=[];

longwavelinear=linspace(min(wavefit(:)),max(wavefit(:)),2338*10);

for o=1:10
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

header=fitsheader(targetfile);
header.IMAGETYP='SPECTRUM';
header.CRVAL1=longwavelinear_clipped(end/2);
header.CRPIX1=length(longwavelinear_clipped)/2;
header.CTYPE1='Wavelength';
header.CUNIT1='nm';
header.CDELT1=mean(diff(longwavelinear_clipped));
headercell=fitstructure2cell(header);
error(' ')
fitswrite(squeeze(finalspeclong_clipped'),[targetfile(1:end-16) '-IndivCalSpec.fit'],headercell(8:end,:))
fitswrite(squeeze(finalspecVarlong_clipped'),[targetfile(1:end-16) '-IndivCalSpec.fit'],'writemode','append')

fitswrite(squeeze(sum(finalspeclong_clipped,2)'),[targetfile(1:end-16) '-CombCalSpec.fit'],headercell(8:end,:))
fitswrite(squeeze(sum(finalspecVarlong_clipped,2)'),[targetfile(1:end-16) '-CombCalSpec.fit'],headercell(8:end,:))
