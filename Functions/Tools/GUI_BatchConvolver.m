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
global selectedMeas selectedFlt outdir %flags
selectedMeas = 0;
selectedFlt = 0;
outdir = 0;
handles.output = hObject;

guidata(hObject, handles);

function varargout = GUI_BatchConvolver_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


% --- Executes on button press in Btn_FilesConvolve.
function Btn_FilesConvolve_Callback(hObject, eventdata, handles)
global selectedMeas
[filename_measurements, pathname_measurements] = uigetfile({'*.wav','WAV-files (*.wav)'},'Pick measurement files','MultiSelect', 'on');

% Check if at least one measurement file was chosen
if iscell(filename_measurements) || ischar(filename_measurements)
    handles.filename_meas = filename_measurements;
    handles.pathname_meas = pathname_measurements;
    selectedMeas = 1;
    guidata(hObject, handles);
end


function Btn_Filter_Callback(hObject, eventdata, handles)
global selectedFlt
[filename_filter, pathname_filter] = uigetfile({'*.wav','WAV-files (*.wav)'},'Pick filter files');

if filename_filter % Same as above...
    handles.filename_filt = filename_filter;
    handles.pathname_filt = pathname_filter;
    selectedFlt = 1;
    guidata(hObject, handles);
end


function Btn_Output_Callback(hObject, eventdata, handles)
global outdir
handles.folder_output = uigetdir();

if handles.folder_output % check if the output folder was chosen
    outdir = 1;
    guidata(hObject, handles);
end


function Btn_Calc_Callback(hObject, eventdata, handles)
global selectedMeas selectedFlt outdir

if (selectedMeas && selectedFlt && outdir)
    [inversefilt, fsIF] = audioread(strcat(handles.pathname_filt, '\', ...
        handles.filename_filt));
    
    convolversystem = dsp.Convolver('Method','Frequency Domain');
    
    if iscell(handles.filename_meas)% check if multiple files where chosen
        
        h = waitbar(0, 'Please wait. This may take several minutes');
        for i = 1:length(handles.filename_meas)
            waitbar(i / length(handles.filename_meas));
            
            [xi, fs] = audioread(strcat( handles.pathname_meas, '\', ...
                char(handles.filename_meas(1, i))));
            
            if fs == fsIF
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
            else
                raiseFsMismatch(char(handles.filename_meas(1, i)), fs, handles.filename_filt, fsIF);
                break
            end
        end
        close(h)
    else % only one file was chosen
        % There is no need for waitbar for only one file... is it?
        
        [x, fs] = audioread(strcat( handles.pathname_meas, '\', ...
                char(handles.filename_meas)));
        
        if fs == fsIF
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
        else
            raiseFsMismatch(char(handles.filename_meas), fs, handles.filename_filt, fsIF);
        end
    end
end


function raiseFsMismatch(file1, fs1, file2, fs2)
    uiwait(msgbox({...
        [ file1 ' sampling frequency is ' num2str(fs1) ' Hz while ' ] 
        [ file2 ' sampling frequency is ' num2str(fs2) ' Hz.' ] ...
    }, 'Sampling frequency mismatch', 'error'));
