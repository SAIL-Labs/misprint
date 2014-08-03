function varargout = wavecalGUI(varargin)
% WAVECALGUI MATLAB code for wavecalGUI.fig
%      WAVECALGUI, by itself, creates a new WAVECALGUI or raises the existing
%      singleton*.
%
%      H = WAVECALGUI returns the handle to a new WAVECALGUI or the handle to
%      the existing singleton*.
%
%      WAVECALGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WAVECALGUI.M with the given input arguments.
%
%      WAVECALGUI('Property','Value',...) creates a new WAVECALGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before wavecalGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to wavecalGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help wavecalGUI

% Last Modified by GUIDE v2.5 17-Jul-2014 17:52:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @wavecalGUI_OpeningFcn, ...
    'gui_OutputFcn',  @wavecalGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before wavecalGUI is made visible.
function wavecalGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wavecalGUI (see VARARGIN)

% Choose default command line output for wavecalGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

order=70:90;
for i=1:length(order)
    contents{i}=num2str(order(i));
end
set(handles.echelleOrder,'String',contents)

guidata(hObject, handles);

% UIWAIT makes wavecalGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = wavecalGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function wave_start_Callback(hObject, eventdata, handles)
% hObject    handle to wave_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wave_start as text
%        str2double(get(hObject,'String')) returns contents of wave_start as a double
try
xlim(handles.refSpectraAxes,[str2double(get(handles.wave_start,'String')) str2double(get(handles.wave_end,'String'))])
catch err
end
% --- Executes during object creation, after setting all properties.
function wave_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wave_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wave_end_Callback(hObject, eventdata, handles)
% hObject    handle to wave_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wave_end as text
%        str2double(get(hObject,'String')) returns contents of wave_end as a double
try
xlim(handles.refSpectraAxes,[str2double(get(handles.wave_start,'String')) str2double(get(handles.wave_end,'String'))])
ylim(handles.refSpectraAxes,[0 max(handles.reference(handles.referenceWave>str2double(get(handles.wave_start,'String'))  & handles.referenceWave<str2double(get(handles.wave_end,'String'))))])
catch err
end

% --- Executes during object creation, after setting all properties.
function wave_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wave_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in orderMenu.
function orderMenu_Callback(hObject, eventdata, handles)
% hObject    handle to orderMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns orderMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from orderMenu
plot(handles.finalSpec(10,:,get(handles.orderMenu,'Value'))/max(handles.finalSpec(10,:,get(handles.orderMenu,'Value'))),'Parent',handles.dataSpectraAxes)
ylim(handles.dataSpectraAxes,[0 1])
xlim(handles.dataSpectraAxes,[1 size(handles.finalSpec,2)])
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function orderMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to orderMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filename_Callback(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename as text
%        str2double(get(hObject,'String')) returns contents of filename as a double


% --- Executes during object creation, after setting all properties.
function filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadData.
function loadData_Callback(hObject, eventdata, handles)
% hObject    handle to loadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



%matdata=load('sun_xcorr_echele_spec_combine.mat','flatflatflatFinalSpec');
%handles.finalSpec=matdata.flatflatflatFinalSpec; clear matdata

if 1
    matdata=load('xcorr_echele_spec_combine.mat','finalSpec');
    handles.finalSpec=matdata.finalSpec; clear matdata
else
    temp=fitsread('sky_cloud-flattend-1D-spectra.fit');
    handles.finalSpec=squeeze(temp(1,:,:));
end
[newSunFlux, newSunWave] = getsunspec(str2double(get(handles.wave_start,'String'))-1,str2double(get(handles.wave_end,'String'))+1,0.03);

if 1
    [TelFlux, TelWave] = getTelluricSpec(str2double(get(handles.wave_start,'String'))-1, str2double(get(handles.wave_end,'String'))+1, 0.03);
    newTelFlux=4.^interp1(TelWave-0.1809,TelFlux,newSunWave);
    newTelFlux=newTelFlux/max(newTelFlux);
    
    newSunFlux=newSunFlux.*newTelFlux;
    
    newSunFlux=newSunFlux/max(newSunFlux);

    clipEdges=newSunWave > str2double(get(handles.wave_start,'String')) & newSunWave < str2double(get(handles.wave_end,'String'));
    newSunFlux=newSunFlux(clipEdges);
    newTelFlux=newTelFlux(clipEdges);
    newSunWave=newSunWave(clipEdges);
    min(newSunWave)
    max(newSunWave)
    handles.newSunFlux=newSunFlux;
    handles.newTelFlux=newTelFlux;
    handles.newSunWave=newSunWave;
    %    roughTelSpecCounterpart=newTelFlux(round([1:length(bestspec')]+info.BestShift));
end

order=get(handles.orderMenu,'Value');


divs=0.685-0.05:0.001:0.68+0.05; %0.042:0.0001:0.048;   %
for i=1:length(divs)
    
    %div=0.02;
    newspec=interp1(1:size(handles.finalSpec,1),handles.finalSpec(:,order),1:divs(i):size(handles.finalSpec,1));
    
    
    plot(newSunWave,newSunFlux,newSunWave,newTelFlux,'r--','Parent',handles.refSpectraAxes)
    axis(handles.refSpectraAxes,[min(newSunWave) max(newSunWave) 0.4 1.1])
    
    plot(1:length(newspec),newspec,'Parent',handles.dataSpectraAxes)
    axis(handles.dataSpectraAxes,[1 length(newspec) 0.4 1.1])
    
    %[c,lags]=xcorr(newflux,newspec,);
    
    [lags,C,info]=xcorr_fft(newSunFlux,newspec','mean');
    
    shifts(i)=info.Shift;
    corrs(i)=info.Corr;
    bestshifts(i)=info.BestShift;
    bestcorrs(i)=info.BestCorr;
    %Info(k)={Infos};
    %plot(lags,C)
    %pause(0.05)
end

%plot(1:length(newspec),newspec+0.1,1:length(newSunFlux),newSunFlux)

[~,ind]=max(bestcorrs);
ind
divs(ind)
shifts(ind)

bestXinterp=1:divs(ind):size(handles.finalSpec,1);
handles.bestXinterp=bestXinterp;
bestspec=interp1([1:size(handles.finalSpec,1)],handles.finalSpec(:,order),bestXinterp);
handles.bestspec=bestspec;

goodwave=round([1:length(bestspec')]+bestshifts(ind));

%goodwave(goodwave>)

bestspec(goodwave<=1 | goodwave>length(newSunWave))=0;
%goodwave(goodwave<=1 | goodwave>length(newSunWave))=0;

roughWave=newSunWave(goodwave);

%min(roughWave)
%max(roughWave)

%max(roughWave)-min(roughWave);

%ax(1)=subplot(3,1,1);


%roughSunSpecCounterpart=newSunFlux(round([1:length(bestspec')]+info.BestShift));


plot(newSunWave,newSunFlux,newSunWave,newTelFlux,'r--','Parent',handles.refSpectraAxes)
axis(handles.refSpectraAxes,[min(roughWave)-0.5 max(roughWave)+0.5 0.4 1.1])
%hold(handles.handles.refSpectraAxes,'all')
grid(handles.refSpectraAxes,'on')

plot(linspace(1,2370,length(roughWave)),bestspec,'Parent',handles.dataSpectraAxes)
axis(handles.dataSpectraAxes,[1 2370 0.4 1.1])
%hold(handles.handles.dataSpectraAxes,'all')
grid(handles.dataSpectraAxes,'on')

plot(roughWave,bestspec-0.1,newSunWave,newSunFlux,newSunWave,newTelFlux,'--','Parent',handles.overlay)
axis(handles.overlay,[min(roughWave)-0.5 max(roughWave)+0.5 0.4 1.1])
%hold(handles.dataSpectraAxes,'all')
grid(handles.overlay,'on')

guidata(hObject, handles);

% --- Executes on button press in fitWave.
function fitWave_Callback(hObject, eventdata, handles)
% hObject    handle to fitWave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for i=1:length(handles.dataLineHandles)
    temp=handles.dataLineHandles(i).getPosition();
    xData(i)=temp(1);
    temp=handles.refLineHandles(i).getPosition();
    xRef(i)=temp(1);
end

xData=sort(xData);
xRef=sort(xRef);

%% find centroids of xref
dataspectra=handles.finalSpec(10,:,get(handles.orderMenu,'Value'))/max(handles.finalSpec(10,:,get(handles.orderMenu,'Value')));

for i=1:length(xData)
    x=round(xData(i))-10:round(xData(i))+10;
    xzoomprofile=dataspectra(x);
    [cf_, xdatafit(i), c, d, a, gof] = fit_gauss(x',xzoomprofile');
end

%%
p=polyfit((xdatafit),(xRef),3)

calwave=polyval(p,1:size(handles.finalSpec,2)); %linspace(1,2370,length(handles.bestXinterp))

if get(handles.saveSpec,'Value')
    fitswrite([handles.bestspec;calwave],['sky_cloud-flattend-1D-spectra-wavecal-order' num2str(get(handles.orderMenu,'Value')) '.fit'])
end

plot(calwave,handles.finalSpec(10,:,get(handles.orderMenu,'Value'))/max(handles.finalSpec(10,:,get(handles.orderMenu,'Value'))),...
    handles.referenceWave,handles.reference,'--','Parent',handles.overlay)

%
%    handles.newSunWave,handles.newSunFlux,'--',...
%    handles.newSunWave,handles.newTelFlux,'--r'

grid(handles.overlay,'on')
axis(handles.overlay,[min(calwave) max(calwave) -0.1 1.1])
filename=get(handles.filename,'String');
save([filename(1:end-5) '-ref-points-order' num2str(get(handles.orderMenu,'Value')) '.mat'],'xData','xRef','p')
guidata(hObject, handles);
%save('temp.mat')

% --- Executes on button press in selectRefPeaks.
function selectRefPeaks_Callback(hObject, eventdata, handles)
% hObject    handle to selectRefPeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%try
    [x,~]=getpts(handles.refSpectraAxes);
    
    contents = cellstr(get(handles.refspecSelection,'String'));
    if strcmpi(contents{get(handles.refspecSelection,'Value')},'ThAr')
        [lines]=findThArline_IRAF(x*10,1)/10;
    elseif strcmpi(contents{get(handles.refspecSelection,'Value')},'CuAr')
        [lines]=findCuAr(x*10,1)/10;
    end
    
    isinlinelist=~isnan(lines);
    
    handles=makeImlinesForEachPoint(handles,lines(isinlinelist),'refLineHandles');
    
    set(handles.wavesTable,'Data',getPosFrom(handles.refLineHandles));
%catch err
%    rethrow(err)
%end
guidata(hObject, handles);

% --- Executes on button press in selectDataPeaks.
function selectDataPeaks_Callback(hObject, eventdata, handles)
% hObject    handle to selectDataPeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[x,~]=getpts(handles.dataSpectraAxes);

handles=makeImlinesForEachPoint(handles,x,'dataLineHandles');

%set(handles.wavesTable,'Data',getPosFrom(handles.dataLineHandles))

guidata(hObject, handles);

function handles=makeImlinesForEachPoint(handles,x,handleNameStr)
switch handleNameStr
    case 'refLineHandles'
        for i=1:length(x)
            
            if isfield(handles,'refLineHandles') & isvalid(handles.refLineHandles)
                handles.refLineHandles(end+1)=imline(handles.refSpectraAxes,[x(i) x(i)],[0 2]);
            else
                handles.refLineHandles(1)=imline(handles.refSpectraAxes,[x(i) x(i)],[0 2]);
                handles.refLineHandles(1)
            end
            addNewPositionCallback(handles.refLineHandles(end),@(pos) updatePosition(pos, handles));
            setColor(handles.refLineHandles(end),'g');
        end
    case 'dataLineHandles'
        for i=1:length(x)
            if isfield(handles,'dataLineHandles') & isvalid(handles.dataLineHandles)
                handles.dataLineHandles(end+1)=imline(handles.dataSpectraAxes,[x(i) x(i)],[0 2]);
            else
                handles.dataLineHandles(1)=imline(handles.dataSpectraAxes,[x(i) x(i)],[0 2]);
            end
            %addNewPositionCallback(handles.refLineHandles(end),@(pos) updatePosition(pos, hObject, handles))
        end
end


% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'dataLineHandles')
    for i=1:length(handles.dataLineHandles)
        try
            handles.dataLineHandles(i).delete();
        catch err
            err
        end
    end
    handles.dataLineHandles(arrayfun(@(s) ~isvalid(s),handles.dataLineHandles))=[];
    handles=rmfield(handles, 'dataLineHandles');
end

if isfield(handles,'refLineHandles')
    for i=1:length(handles.refLineHandles)
        try
            handles.refLineHandles(i).delete();
        catch err
            err
        end
    end
    handles.refLineHandles(arrayfun(@(s) ~isvalid(s),handles.refLineHandles))=[];
    handles=rmfield(handles, 'refLineHandles');
end

set(handles.wavesTable,'Data',[])

%%catch
%%end
guidata(hObject, handles);

function xData=getPosFrom(LineHandles)
for i=1:length(LineHandles)
    temp=LineHandles(i).getPosition();
    xData(i)=temp(1);
end
xData=sort(xData');


function updatePosition(pos, handles)
handles=guidata(handles.output); % to get current verions
set(handles.wavesTable,'Data',getPosFrom(handles.refLineHandles))
guidata(handles.output, handles);



% --- Executes during object creation, after setting all properties.
function wavesTable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wavesTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'ColumnName','Reference Wave','Data',{blanks(0);blanks(0);blanks(0)});


% --- Executes on button press in loadSavedRefPoints.
function loadSavedRefPoints_Callback(hObject, eventdata, handles)
% hObject    handle to loadSavedRefPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename=get(handles.filename,'String');

savedata=load([filename(1:end-5) '-ref-points-order' num2str(get(handles.orderMenu,'Value')) '.mat']);

handles=makeImlinesForEachPoint(handles,savedata.xData,'dataLineHandles');
handles=makeImlinesForEachPoint(handles,savedata.xRef,'refLineHandles');
guidata(hObject, handles);

% --- Executes on button press in saveSpec.
function saveSpec_Callback(hObject, eventdata, handles)
% hObject    handle to saveSpec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in refspecSelection.
function refspecSelection_Callback(hObject, eventdata, handles)
% hObject    handle to refspecSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns refspecSelection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from refspecSelection
contents = cellstr(get(hObject,'String'));

if strcmpi(contents{get(hObject,'Value')},'ThAr')
    [Wave,Intensity]=ThArSpec(5400,8000); %str2double(get(handles.wave_start,'String'))*10,str2double(get(handles.wave_end,'String'))*10
    handles.reference=Intensity/max(Intensity);
    handles.referenceWave=Wave/10;
elseif strcmpi(contents{get(hObject,'Value')},'CuAr')
    
    [cuarspec, header]=fitsread('cuar.fits');
    wave=header.CRVAL1+header.CDELT1*((1:length(cuarspec))-header.CRPIX1);
    
    handles.reference=cuarspec/max(cuarspec);
    handles.referenceWave=wave/10;
end
plot(handles.referenceWave,handles.reference,'Parent',handles.refSpectraAxes)
xlim(handles.refSpectraAxes,[str2double(get(handles.wave_start,'String')) str2double(get(handles.wave_end,'String'))])
ylim(handles.refSpectraAxes,[-0.1 1.1])
grid(handles.refSpectraAxes,'on')
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function refspecSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to refspecSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in openAndSelectFileButton.
function openAndSelectFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to openAndSelectFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,PathName] = uigetfile({'*.fits';'*.fit';'*.mat'},'File Selector');
if filename
    cd(PathName);
    handles.finalSpec=fitsread(filename);
    set(handles.filename,'String',filename);
    
    plot(handles.finalSpec(10,:,get(handles.orderMenu,'Value'))/max(handles.finalSpec(10,:,get(handles.orderMenu,'Value'))),'Parent',handles.dataSpectraAxes)
    ylim(handles.dataSpectraAxes,[0 1.1])
    xlim(handles.dataSpectraAxes,[1 size(handles.finalSpec,2)])
    grid(handles.dataSpectraAxes,'on')
end
guidata(hObject, handles);
 


% --- Executes on selection change in echelleOrder.
function echelleOrder_Callback(hObject, eventdata, handles)
% hObject    handle to echelleOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns echelleOrder contents as cell array
%        contents{get(hObject,'Value')} returns selected item from echelleOrder
try
contents = cellstr(get(hObject,'String'));
m=str2double(contents{get(hObject,'Value')});
centralWave=2*(1e-3/31.6)*cosd(5)*sind(63.2)./m / 1e-9;

set(handles.wave_start,'String',num2str(centralWave-abs(centralWave/m/2)*1.1))
set(handles.wave_end,'String',num2str(centralWave+abs(centralWave/m/2)*1.1))

xlim(handles.refSpectraAxes,[str2double(get(handles.wave_start,'String')) str2double(get(handles.wave_end,'String'))])
catch err
    err
end


% --- Executes during object creation, after setting all properties.
function echelleOrder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to echelleOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
