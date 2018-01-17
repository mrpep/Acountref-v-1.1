function varargout = GUI_BatchConvolver(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_BatchConvolver_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_BatchConvolver_OutputFcn, ...
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

function GUI_BatchConvolver_OpeningFcn(hObject, eventdata, handles, varargin)
global openedmeas openedfilt outdir
openedmeas = 0;
openedfilt = 0;
outdir = 0;
handles.output = hObject;

guidata(hObject, handles);

function varargout = GUI_BatchConvolver_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


% --- Executes on button press in Btn_FilesConvolve.
function Btn_FilesConvolve_Callback(hObject, eventdata, handles)
global openedmeas
[filename_measurements, pathname_measurements] = uigetfile({'*.wav','WAV-files (*.wav)'},'Pick measurement files','MultiSelect', 'on');

handles.filename_meas = filename_measurements;
handles.pathname_meas = pathname_measurements;

openedmeas = 1;
guidata(hObject, handles);

function Btn_Filter_Callback(hObject, eventdata, handles)
global openedfilt
[filename_filter, pathname_filter] = uigetfile({'*.wav','WAV-files (*.wav)'},'Pick filter files');

handles.filename_filt = filename_filter;
handles.pathname_filt = pathname_filter;

openedfilt = 1;
guidata(hObject, handles);

function Btn_Output_Callback(hObject, eventdata, handles)
global outdir
handles.folder_output = uigetdir();
outdir = 1;
guidata(hObject, handles);

function Btn_Calc_Callback(hObject, eventdata, handles)
global openedmeas openedfilt outdir

if (openedmeas && openedfilt && outdir)
    [inversefilt, fs] = audioread(strcat(handles.pathname_filt, '\', ...
        handles.filename_filt));
    
    convolversystem = dsp.Convolver('Method','Frequency Domain');
    
    if iscell(handles.filename_meas)% check if multiple files where chosen
        
        h = waitbar(0, 'Please wait. This may take several minutes');
        for i = 1:length(handles.filename_meas)
            waitbar(i / length(handles.filename_meas));
            
            [xi, fs] = audioread(strcat( handles.pathname_meas, '\', ...
                char(handles.filename_meas(1, i))));
            
            impulse = step(convolversystem, xi, inversefilt);
            
            if (get(handles.Radio_Normalize, 'Value'))
                impulse = impulse / max(abs(impulse));
            elseif (get(handles.Radio_Preserve, 'Value'))
                if i == 1
                    divisor = max(abs(impulse));
                end
                impulse = 0.8*(impulse/divisor);
            end
            
            audiowrite(strcat(handles.folder_output, '\', ...
                char(handles.filename_meas(1, i))), impulse, fs);
        end
        close(h)
    else % only one file was chosen
        % There is no need for waitbar for only one file... is it?
        
        [x, fs] = audioread(strcat( handles.pathname_meas, '\', ...
                char(handles.filename_meas)));
        
        impulse = step(convolversystem, x, inversefilt);
            
        if (get(handles.Radio_Normalize, 'Value'))
            impulse = impulse / max(abs(impulse));
        elseif (get(handles.Radio_Preserve, 'Value'))
            divisor = max(abs(impulse));
            impulse = 0.8*(impulse/divisor);
        end
        
        audiowrite(strcat(handles.folder_output, '\', ...
                char(handles.filename_meas)), impulse, fs);
        
        uiwait(msgbox([ ...
            'Done. It''s named "Batch Convolver" for a reason, you ' ...
            'know... you can convolve multiple files by holding the ' ...
            'Ctrl key while picking them.' ...
        ], 'Convolution done', 'modal'))
    end
end
