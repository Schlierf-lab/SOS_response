% User interface for interactive cell selection
% --Andreas Hartmann

function varargout = CellSelector(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CellSelector_OpeningFcn, ...
                   'gui_OutputFcn',  @CellSelector_OutputFcn, ...
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

function CellSelector_OpeningFcn(hObject, eventdata, handles, varargin)

global plotHandle;

handles.output = hObject;

guidata(hObject, handles);

set(handles.axesBF,'Box','on');
set(handles.axesBF,'XTickLabel',{});
set(handles.axesBF,'YTickLabel',{});

plotHandle=[];

function varargout = CellSelector_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

function listbox1_Callback(hObject, eventdata, handles)

function listbox1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Untitled_1_Callback(hObject, eventdata, handles)

function FOpen_Callback(hObject, eventdata, handles)

global FILE;
global PATH;
global cellborders;
global cellborders_unCorr;

[FILE,PATH]=uigetfile('*.tif');

% get reader for the nd2 files
imageBF=imread([PATH FILE]);

% load cell borders
cellborders_unCorr=matfile([PATH FILE(1:end-7) '_out.mat'],'Writable',true);
cellborders=matfile([PATH FILE(1:end-7) '_out_corr.mat'],'Writable',true);
loaddata=cellborders.cellList;
loaddata_unCorr=cellborders_unCorr.cellList;

% display brightfield image
axes(handles.axesBF);
imagesc(handles.axesBF,imageBF);
colormap('gray');
hold(handles.axesBF,'on');

% clean meshdata
meshData={[]};
meshData{1}={[]};
meshData_unCorr={[]};
meshData_unCorr{1}={[]};

iterN=1;
for iter=1:cellborders.cellListN

    if ~isempty(loaddata.meshData{1}{iter}.model)
        
        meshData{1}{iterN}=loaddata.meshData{1}{iter};
        meshData_unCorr{1}{iterN}=loaddata_unCorr.meshData{1}{iter};
        iterN=iterN+1;
    end
end

loaddata.meshData=meshData;
loaddata_unCorr.meshData=meshData_unCorr;

cellborders.cellList=loaddata;
cellborders.cellListN=iterN-1;
cellborders_unCorr.cellList=loaddata_unCorr;
cellborders_unCorr.cellListN=iterN-1;

% load data into list
data=cell(cellborders.cellListN,3);

for iter=1:cellborders.cellListN

    str=['Cell_' num2str(iter)];
    data{iter,1}=str;
    
    if isfield(loaddata,'selected1')&&isfield(loaddata,'selected2')
    
        data{iter,2}=loaddata.selected1(iter)==1;
        data{iter,3}=loaddata.selected2(iter)==1;
    else
        data{iter,2}=true;
        data{iter,3}=true;
    end
end

set(handles.uitable,'Data',data);
set(handles.textFile,'String',FILE(1:end-4));

function uitable_CellSelectionCallback(hObject, eventdata, handles)

global plotHandle;
global cellborders_unCorr;

entries=get(handles.uitable,'Data');
load('cmap.txt');

if ~isempty(plotHandle)
    
    delete(plotHandle);
end

vIndices=eventdata.Indices;

if ~isempty(vIndices)

    data=cellborders_unCorr.cellList;
    
    if vIndices(2)==1

        plotHandle=plot(handles.axesBF,data.meshData{1}{vIndices(1)}.model(:,1),data.meshData{1}{vIndices(1)}.model(:,2),'LineWidth',2);
        
        set(plotHandle,'Color',cmap(entries{vIndices(1),2}+entries{vIndices(1),3}*2+4,:));
    end
end

function Save_Callback(hObject, eventdata, handles)

global cellborders;

entries=get(handles.uitable,'Data');

arrSel1=zeros(1,size(entries,1));
arrSel2=zeros(1,size(entries,1));

h=waitbar(0,'Please wait...');

for iter=1:size(entries,1)
    
    arrSel1(iter)=entries{iter,2};
    arrSel2(iter)=entries{iter,3};
    
    waitbar(iter/size(entries,1),h);
end

close(h);

data=cellborders.cellList;
data.selected1=arrSel1;
data.selected2=arrSel2;

cellborders.cellList=data;
