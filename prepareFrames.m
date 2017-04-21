%function prepareFrames
    %preparedframes
    %clc; clear
    fits='fits';
    fits='fit';

    
    %% scan for darks, light exposure times
    files=dir(['*.' fits]);
    allfilenames={files.name}';
    allfilenames=allfilenames(cellfun('isempty', strfind(allfilenames,'reduced')));
    
    if isempty(allfilenames)
        error('No Fits files')
    end
    %%
    [isbias,isdark,ismasterdark,ismasterbias,islight,isflat]=deal(false(1,length(allfilenames)));
    imagetypekeywords={'OBSTYPE','IMAGETYP'};
    
    parfor i=1:length(allfilenames)
        try
            header=fitsheader(allfilenames{i});
            try
                exposure(i)=header.EXPOSURE;
            catch
                exposure(i)=header.EXPOSED;
            end
            
            if headerKeywordValueCheck( header, imagetypekeywords, {'BIAS'})
                isbias(i)=true;
            elseif headerKeywordValueCheck( header, imagetypekeywords , {'Dark Frame','DARK'})
                isdark(i)=true;
            elseif headerKeywordValueCheck( header, imagetypekeywords, {'MASTDARK'})
                ismasterdark(i)=true;
            elseif headerKeywordValueCheck( header, imagetypekeywords, {'MASTBIAS'})
                ismasterbias(i)=true;
            elseif headerKeywordValueCheck( header, imagetypekeywords, {'LIGHT','Light Frame','ARC','OBJECT'})
                islight(i)=true;
            elseif headerKeywordValueCheck( header, imagetypekeywords, {'Flat','Flat Field'})
                isflat(i)=true;
            end
        catch err
            allfilenames{i}
            header
            rethrow(err)
        end
    end
    
    processedfolder='processed';
    if ~exist( processedfolder, 'dir' )
        mkdir( processedfolder );
    end
    
    
    %% bias
    try
        bias=fitsread(['masterbias.' fits]);
    catch err   
        % prepare bias frame
        biasfilenames={allfilenames{isbias}}';
        if ~isempty(biasfilenames)
            currentProcceesedFolder=fullfile(processedfolder,['bias']);
            
            if ~exist( currentProcceesedFolder, 'dir' )
                mkdir( currentProcceesedFolder);
            end
            
            
            header=fitsheader(biasfilenames{1});
            header=setMutipleFITSKeywords( header,imagetypekeywords,'MASTBIAS');
            headercell=fitstructure2cell(header);
            
            clear frames
            for ii=1:length(biasfilenames)
                frames(:,:,ii)=uint16(fitsread(biasfilenames{ii}));
                movefile(biasfilenames{ii},currentProcceesedFolder)
            end
            frame=trimmean(frames,10,3);
            
            fitswrite(single(frame),['masterbias.' fits],headercell)
            
        else
            disp('No master Bias and No Bias Frames to process. Bias set to Zero')
            bias=0
        end
    end
    
    %% dark
    % check all required darks are avaliable
    lightexposures=unique(exposure(exposure~=0));
    darkexposures=unique(exposure(isdark | ismasterdark));
    requiredDarks=setdiff(lightexposures,darkexposures);
    
    assertWarn(~isempty(requiredDarks),'MISPRINT:prepareFrame:darksFramesMissing',['These exposure times have no corresponding darks: ' num2str(requiredDarks)])
    
    % preparedark frame, based on exposure, and filename
    darkfilename={allfilenames{isdark}}';
    if ~isempty(darkfilename)
        % loop through each unique dark
        for i=1:length(darkexposures)
            filenames=allfilenames(exposure==darkexposures(i) & isdark);
            currentProcceesedFolder=fullfile(processedfolder,['darks-' num2str(darkexposures(i))]);
            
            if ~exist( currentProcceesedFolder, 'dir' )
                mkdir( currentProcceesedFolder);
            end
            try
                fitsread([num2str(darkexposures(i)) '-masterdark.' fits]);
            catch err
                
                header=fitsheader(filenames{1});
                header.IMAGETYP='MASTDARK';
                headercell=fitstructure2cell(header);
                
                clear frames
                for ii=1:length(filenames)
                    frames(:,:,ii)=uint16(fitsread(filenames{ii}))-bias;
                    movefile(filenames{ii},currentProcceesedFolder)
                end
                frame=mean(frames,3);
                
                fitswrite(single(frame),[num2str(darkexposures(i)) '-masterdark.' fits],headercell)
            end
        end
    else
        disp('No Darks to process.')
    end
    
    %% load masterdarks
    masterdarkfilename={allfilenames{ismasterdark}}';
    if ~isempty(masterdarkfilename)
        for i=1:length(masterdarkfilename)
            header=fitsheader(masterdarkfilename{i});
            darks.(['D' num2str(header.EXPOSURE*100)])=fitsread(masterdarkfilename{i});
        end
    else
        return
        % normally always want darks
        for i=1:length(lightexposures)
            if lightexposures(i)>=1
                darks.(['D' num2str(lightexposures(i)*100)])=0;
            end
        end
    end
    
    %darks.D100=darks.D60;
    
    %% prepare light frames
    lightfilename={allfilenames{islight}}';
    ligthexposures=exposure(islight);
    
    if ~isempty(lightfilename)
        
        %     for i=1:length(lightfilename)
        %         lightbasefilenames(i)={lightfilename{i}(1:end-8)}';
        %     end
        %
        %     lightbasefilenames=unique(lightbasefilenames);
        
        if ~exist( fullfile(processedfolder,'light'), 'dir' )
            mkdir( fullfile(processedfolder,'light') );
        end
        % loop through each unique dark
        parfor i=1:length(lightfilename)
            try
                fitsread([lightfilename{i}(1:end-4) '-reduced.fit']);
            catch err
                header=fitsheader(lightfilename{i});
                header=setMutipleFITSKeywords( header,imagetypekeywords,'reduced');
                headercell=fitstructure2cell(header);
                try
                    %clear frame
                    frame=fitsread(lightfilename{i}) - darks.(['D' num2str(ligthexposures(i)*100)])-bias;
                    %frame(frame<0)=0;
                    movefile(lightfilename{i},fullfile(processedfolder,'light'))
                    fitswrite(single(frame),[stripextension(lightfilename{i}) '-reduced.' fits],headercell)
                    
                catch err
                    if strcmpi(err.identifier,'MATLAB:nonExistentField')
                        warning([lightfilename{i} ': No dark with correct exposure time (probably)'])
                    else
                        rethrow(err)
                    end
                end
                
                
            end
        end
    else
        disp('No Lights to process.')
    end
    
    
    %%proccess flat
    
    flatfilename={allfilenames{isflat}}';
    flatexposures=exposure(isflat);
    if ~isempty(flatfilename)
        
        %     for i=1:length(lightfilename)
        %         lightbasefilenames(i)={lightfilename{i}(1:end-8)}';
        %     end
        %
        %     lightbasefilenames=unique(lightbasefilenames);
        
        if ~exist( fullfile(processedfolder,'flat'), 'dir' )
            mkdir( fullfile(processedfolder,'flat') );
        end
        % loop through each unique dark
        for i=1:length(flatfilename)
            try
                fitsread([flatfilename{i}(1:end-4) '-reduced.' fits]);
            catch err
                header=fitsheader(flatfilename{i});
                header=setMutipleFITSKeywords(header,imagetypekeywords,'reduced');
                headercell=fitstructure2cell(header);
                try
                    clear frame
                    frame=fitsread(flatfilename{i}) - darks.(['D' num2str(flatexposures(i)*100)])-bias;
                    %frame(frame<0)=0;
                    movefile(flatfilename{i},fullfile(processedfolder,'flat'))
                    fitswrite(single(frame),[stripextension(flatfilename{i}) '-reduced.' fits],headercell)
                    
                catch err
                    if strcmpi(err.identifier,'MATLAB:nonExistentField')
                        warning([flatfilename{i} ': No dark with correct exposure time (probably)'])
                    else
                        rethrow(err)
                    end
                end
            end
        end
    else
        disp('No Flats to process.')
    end
    
    return
    %%
    
    for i=1:20
        %for j=[60]
        %movefile(['dark-' num2str(i,'%.3d') '-' num2str(j) '.fit'],['D-' num2str(j) '-' num2str(i,'%.3d') '-dark.fit'])
        %movefile(['dark-' num2str(j) '-' num2str(i,'%.3d') '.fit'],['dark-' num2str(j) '-' num2str(i,'%.3d') '-dark.fit'])
        %movefile(['lamp-' num2str(i,'%.3d') '-60.fit'],['lamp-' num2str(i,'%.3d') '.fit'])
        %end
        movefile(['dark-' num2str(i,'%.3d') '.fit'],['dark-' num2str(i,'%.3d') '-dark.fit'])
    end
    
    
    
    return
    %%
    
    for i=1:5
        lightdark(:,:,i)=fitsread(['dark-' num2str(i,'%.3d') 'ld-reduced.fit']);
    end
    
    lightdark=median(lightdark,3);
    
    %%
    for i=1:10
        %     [newsun header]=fitsread(['flat-' num2str(i,'%.3d') '-reduced.fit']);
        %
        %     header.IMAGETYP='reduced';
        %     headercell=fitstructure2cell(header);
        %
        %     newsun=newsun-lightdark;
        %     newsun=newsun+abs(min2(newsun));
        %
        %     fitswrite(newsun,['flat-' num2str(i,'%.3d') '-reduced-final.fit'],headercell(8:end,:))
        %fitsAddHeaderKeyword(['sun-' num2str(i,'%.3d') '-reduced-final.fit'],'DISPAXIS',1,' ')
        %fitsAddHeaderKeyword(['dark-' num2str(i,'%.3d') '-ldark.fit'],'IMAGETYP','DARK',' ')
    end
    return
    %%
    for i=1:10
        movefile(['flat-' num2str(i,'%.3d') 'dark.fit'],['flat-' num2str(i,'%.3d') '-dark.fit'])
    end
    
    %%
    files=dir('*.fit')
    filenames={files.name}';
    filedate=datenum({files.date});
    [~,idx]=sort(filedate);
    filenames=filenames(idx);
    
    dirFilenames('HD153135*.fit')
    
    for i=1:length(filenames)
        movefile(filenames{i},['oops-' num2str(i,'%.3d') '.fit'])
    end
    
    
    for i=1:51
        imdata(:,:,i)=fitsread(['hene-' num2str(i,'%.3d') '-reduced.fit']);
    end
    %%
    
    
    imdata2=bsxfun(@minus,imdata,median(imdata,3));
    imdata2(imdata2<0)=0;
    %%
    imagesc(sum(imdata2,3))
    %%
    fitswrite(max(imdata2,[],3),'maxCombhene.fit')
    fitswrite(sum(imdata2,3),'sumCombhene.fit')
    
    
    
    %%
    
    filenames=dirFilenames('*-dark30.fit')
    
    for i=1:length(filenames)
        %disp(['oops-' num2str(i,'%.3d') '.fit'])
        movefile(filenames{i},['oops-' num2str(i,'%.3d') '-dark.fit'])
    end
    
