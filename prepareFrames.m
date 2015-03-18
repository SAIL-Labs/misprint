function prepareFrames
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
    fitswrite(single(bias),'master-bias.fit')
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
        if strcmpi(header.IMAGETYP,'DARK') || strcmpi(header.IMAGETYP,'Dark Frame')
            isdark(i)=true;
            %ismaster(i)=true;
        elseif strcmpi(header.IMAGETYP,'MASTDARK')
            ismasterdark(i)=true;
        elseif strcmpi(header.IMAGETYP,'LIGHT') || strcmpi(header.IMAGETYP,'LIGHT Frame')
            islight(i)=true;
        end
    catch
        allfilenames{i}
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
            
            fitswrite(single(frame),[basefilename '-masterdark.fit'],headercell(8:end,:))
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

%darks.D100=darks.D60;


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
                frame=fitsread(lightfilename{i}) - darks.(['D' num2str(header.EXPOSURE*100)]);
                %frame(frame<0)=0;
                movefile(lightfilename{i},fullfile(processedfolder,lightfilename{i}(1:end-8)))
                fitswrite(single(frame),[lightfilename{i}(1:end-4) '-reduced.fit'],headercell(8:end,:))
            
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
%%
% for i=1:length(filenames)
%     movefile(filenames{i},['oops-' num2str(i,'%.3d') '.fit'])
% end
% 

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







