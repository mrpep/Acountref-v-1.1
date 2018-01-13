function varargout = GUI_BatchCutter(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_BatchCutter_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_BatchCutter_OutputFcn, ...
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

function GUI_BatchCutter_OpeningFcn(hObject, eventdata, handles, varargin)

handles.NSegments = 0;
handles.Segments = cell(1,3);
handles.output = hObject;
guidata(hObject, handles);

function varargout = GUI_BatchCutter_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

function Btn_Open_Callback(hObject, eventdata, handles)

[filename_audios, pathname_audios] = uigetfile({'*.wav','WAV-files (*.wav)'},'Pick audio files','MultiSelect', 'on');

handles.filename_audios = filename_audios;
handles.pathname_audios = pathname_audios;

[audiosample,fs] = audioread(strcat(pathname_audios,'\',char(filename_audios(1,1))));
audiosample = audiosample/max(abs(audiosample));
t = (0:length(audiosample)-1)/fs;
plot(handles.Plot_Audio,t,audiosample)

guidata(hObject, handles);

function Btn_Add_Callback(hObject, eventdata, handles)

NSegments = handles.NSegments;
NSegments = NSegments + 1;
Segments = handles.Segments;
prompt = {'Enter segment name:'};
dlg_title = 'New Segment';
num_lines = 1;
defaultans = {num2str(handles.NSegments)};
segmentname = inputdlg(prompt,dlg_title,num_lines,defaultans);
%Enable mouse selection of begin and end of segment
[xbegin,~] = ginput(1);
hold(handles.Plot_Audio,'on');
stem(handles.Plot_Audio,xbegin,1,'Marker','none','Color','g');
[xend,~] = ginput(1);
hold(handles.Plot_Audio,'on');
stem(handles.Plot_Audio,xend,1,'Marker','none','Color','r');

if (xbegin<xend)
    Segments{NSegments,1} = segmentname;
    Segments{NSegments,2} = xbegin;
    Segments{NSegments,3} = xend;
else
    Segments{NSegments,1} = segmentname;
    Segments{NSegments,2} = xend;
    Segments{NSegments,3} = xbegin;
end
handles.Segments = Segments;
handles.NSegments = NSegments;

guidata(hObject, handles);

function Btn_Save_Callback(hObject, eventdata, handles)

outdir = uigetdir();
Segments = handles.Segments;
filenames_audio = handles.filename_audios;
[a,~] = size(Segments);
h = waitbar(0,'Please wait. This may take some seconds');
for i = 1:length(filenames_audio)
    waitbar(i/length(filenames_audio))
    [audiosample,fs] = audioread(strcat(handles.pathname_audios,'\',char(filenames_audio(1,i))));
    for j = 1:a
        kbegin = round(fs*cell2mat(Segments(j,2)));
        kend = round(fs*cell2mat(Segments(j,3)));
        if (kbegin<1)
            kbegin = 1;
        end
        if (kend>length(audiosample))
            kend = length(audiosample);
        end        
        audiosegment = audiosample(kbegin:kend);
        audiowrite(strcat(outdir,'\',char(Segments{j,1}(1,1)),'_',char(filenames_audio(1,i))),audiosegment,fs);        
    end
end
close(h)
