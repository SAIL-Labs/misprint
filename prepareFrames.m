%preparedframes
clc; clear all

processedfolder='processed';
if ~exist( processedfolder, 'dir' )
    mkdir( processedfolder );
end

%% bias
basefilename='bias';
files=dir([basefilename '-*.fit']);
filenames={files.name}';

if ~exist( fullfile(processedfolder,basefilename), 'dir' )
    mkdir( fullfile(processedfolder,basefilename) );
end
try
    bias=fitsread('master-bias.fit');
catch err
    try 
    for i=1:length(filenames)
        biasdata(:,:,i)=fitsread(filenames{i});
        movefile(filenames{i},fullfile(processedfolder,basefilename))
    end
    bias=median(biasdata,3);
    fitswrite(bias,'master-bias.fit')
    catch err
        bias=0;
    end
end


%% scan for darks, light exposure times
files=dir('*.fit');
allfilenames={files.name}';
%header=[];
isdark=logical(zeros(1,length(allfilenames)));
ismasterdark=logical(zeros(1,length(allfilenames)));
islight=logical(zeros(1,length(allfilenames)));
for i=1:length(allfilenames)
    try
        header=fitsheader(allfilenames{i});
        exposure(i)=header.EXPOSURE;
        if strcmpi(header.IMAGETYP,'DARK')
            isdark(i)=true;
            %ismaster(i)=true;
        elseif strcmpi(header.IMAGETYP,'MASTDARK')
            ismasterdark(i)=true;
        elseif strcmpi(header.IMAGETYP,'LIGHT')
            islight(i)=true;
        end
    catch
        header
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
    for i=1:length(darkfilename)
        darkbasefilenames(i)={darkfilename{i}(1:end-13)}';
    end
    
    darkbasefilenames=unique(darkbasefilenames);
    
    % loop through each unique dark
    for i=1:length(darkbasefilenames)
        basefilename=darkbasefilenames{i};
        files=dir([basefilename '-*-dark.fit']);
        filenames={files.name}';
        
        if ~exist( fullfile(processedfolder,basefilename), 'dir' )
            mkdir( fullfile(processedfolder,basefilename) );
        end
        try
            fitsread([darkbasefilenames '-masterdark.fit']);
        catch err
            
            header=fitsheader(filenames{1});
            header.IMAGETYP='MASTDARK';
            headercell=fitstructure2cell(header);
            
            clear frames
            for i=1:length(filenames)
                frames(:,:,i)=fitsread(filenames{i})-bias;
                movefile(filenames{i},fullfile(processedfolder,basefilename))
            end
            frame=median(frames,3);
            
            fitswrite(frame,[basefilename '-masterdark.fit'],headercell(8:end,:))
        end
    end
else
    disp('No Darks to process.')
end

%% load masterdarks
darkfilename={allfilenames{ismasterdark}}';
for i=1:length(darkfilename)
    header=fitsheader(darkfilename{i});
    darks.(['D' num2str(header.EXPOSURE*100)])=fitsread(darkfilename{i});
end

%% prepare light frames
lightfilename={allfilenames{islight}}';

if ~isempty(lightfilename)
    
    %     for i=1:length(lightfilename)
    %         lightbasefilenames(i)={lightfilename{i}(1:end-8)}';
    %     end
    %
    %     lightbasefilenames=unique(lightbasefilenames);
    
    
    % loop through each unique dark
    for i=1:length(lightfilename)
        
        
        if ~exist( fullfile(processedfolder,lightfilename{i}(1:end-8)), 'dir' )
            mkdir( fullfile(processedfolder,lightfilename{i}(1:end-8)) );
        end
        try
            fitsread([lightfilename{i}(1:end-4) '-reduced.fit']);
        catch err
            
            header=fitsheader(lightfilename{i});
            header.IMAGETYP='reduced';
            headercell=fitstructure2cell(header);
            try
                clear frame
                frame=fitsread(lightfilename{i}) - bias - darks.(['D' num2str(header.EXPOSURE*100)]);
                frame(frame<0)=0;
                movefile(lightfilename{i},fullfile(processedfolder,lightfilename{i}(1:end-8)))
                fitswrite(frame,[lightfilename{i}(1:end-4) '-reduced.fit'],headercell(8:end,:))
            
            catch err
                err
                warning([lightfilename{i} ': No dark with correct exposure time (probably)'])
            end
            
            
        end
    end
else
    disp('No Lights to process.')
end

return
%%

for i=1:20
    %for j=[60]
        %movefile(['dark-' num2str(i,'%.3d') '-' num2str(j) '.fit'],['D-' num2str(j) '-' num2str(i,'%.3d') '-dark.fit'])
        %movefile(['dark-' num2str(j) '-' num2str(i,'%.3d') '.fit'],['dark-' num2str(j) '-' num2str(i,'%.3d') '-dark.fit'])
        %movefile(['lamp-' num2str(i,'%.3d') '-60.fit'],['lamp-' num2str(i,'%.3d') '.fit'])
    %end
    movefile(['dark-' num2str(i,'%.3d') 'dark.fit'],['dark-' num2str(i,'%.3d') '-dark.fit'])
end











