%Positions Cell Array (i,j,k):
% i= pos
% j =
    % 1 = Name
    % 2 = X
    % 3 = Y
    % 4 = Mic Number
    % 5 = Type: 1 = Monoaural 2 = Binaural 3 = Soundfield
    % 6 = IR File
    % 7 = WN File
    % 8-11 = Anechoic 1-4 File
    % 12 = Calibration Factor (94 dB)
% k =  Fuente

% SF:

% Archivo Stereo, L = Figura 8, R = Omni

%Cambiar enfoque de IR integradas subjetivamente -> calcularlas una vez de
%entrada y listo.

function varargout = GUI_Main(varargin)

    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @GUI_Main_OpeningFcn, ...
                       'gui_OutputFcn',  @GUI_Main_OutputFcn, ...
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

function GUI_Main_OpeningFcn(hObject, eventdata, handles, varargin)

    %Se borra basura previa al abrir y se inicializa el indice de posiciones:
    clearvars -global
    global npositions Positions counter NSources
    addpath(strcat(pwd,'\Functions'));
    addpath(strcat(pwd,'\Functions\Parameters'));
    addpath(strcat(pwd,'\Functions\Tasks'));
    addpath(strcat(pwd,'\Functions\Tools'));
    addpath(strcat(pwd,'\Functions\Utils'));
    npositions = 0; %numero de posiciones
    counter = 0; %contador para el nombre de posicion
    NSources = 2; %Cambiar si se quieren más posiciones de fuente
    for i=1:NSources       
        SourcesStr(i) = cellstr(num2str(i));
    end
    SourcesStr(NSources + 1) = cellstr('Mean');
    set(handles.Combo_SourceMap,'String',SourcesStr);
    set(handles.Combo_SourceSet,'String',SourcesStr(1:NSources));
       
    handles.ScaleDefined = 0;
    handles.output = hObject;
    guidata(hObject,handles);
    warning('off');
    
    uiwait(msgbox('First click on the button Load Floor Plan and select an image of the room you measured (No smoke allowed)','Welcome to Acuntref','modal'))

    function varargout = GUI_Main_OutputFcn(hObject, eventdata, handles) 

    varargout{1} = handles.output;

function Combo_Parameter_Callback(hObject, eventdata, handles)
    
    sel = get(handles.Combo_Parameter,'Value');
    if(sel<27)
        col2use = handles.colstr;
    else
        col2use = handles.col2str;
    end
    set(handles.Combo_Band,'String',col2use);
    
function Combo_Parameter_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function Combo_Band_Callback(hObject, eventdata, handles)

function Combo_Band_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function Combo_Colormap_Callback(hObject, eventdata, handles)

function Combo_Colormap_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function List_Positions_Callback(hObject, eventdata, handles)
    %Al hacer clic en la lista se muestra la info de esa posición.
    global Positions ResultsObtained
    source = get(handles.Combo_SourceSet,'Value');
    selpos = get(handles.List_Positions,'Value');
    set(handles.Edit_XPix,'String',num2str(round(cell2mat(Positions(selpos,2,1)))));
    set(handles.Edit_YPix,'String',num2str(round(cell2mat(Positions(selpos,3,1)))));
    PlotPositions(handles,hObject,selpos);
    if (handles.ScaleDefined==1)
        Xpix = str2num(get(handles.Edit_XPix,'String'));
        Ypix = str2num(get(handles.Edit_YPix,'String'));
        XMeters = Xpix*handles.scale;
        YMeters = (handles.imagesize(1)-Ypix)*handles.scale;
        set(handles.Edit_XMeters,'String',num2str(XMeters));
        set(handles.Edit_YMeters,'String',num2str(YMeters));  
    end
    strlist = get(handles.List_Positions,'String');
    set(handles.Edit_PositionName,'String',strlist(selpos));
    set(handles.Edit_MicNumber,'String',num2str(round(cell2mat(Positions(selpos,4,1)))));
    if(Positions{selpos,5,1}>0)
        set(handles.Combo_PosType,'Value',Positions{selpos,5,1});
    else
        set(handles.Combo_PosType,'Value',1);
    end
    set(handles.Edit_IRFile,'String',Positions{selpos,6,source});
    set(handles.Edit_WNFile,'String',Positions{selpos,7,source});
    set(handles.Edit_AN1File,'String',Positions{selpos,8,source});
    set(handles.Edit_AN2File,'String',Positions{selpos,9,source});
    set(handles.Edit_AN3File,'String',Positions{selpos,10,source});
    set(handles.Edit_AN4File,'String',Positions{selpos,11,source});
    cla(handles.Plot_IR)
    
    %Si ya se cargo IR se grafica:
    if (isempty(Positions{selpos,6,source})==0)
        [selIR,fs] = audioread(char(Positions(selpos,6,source)));
        timeIR = (0:length(selIR)-1)/fs;
        plot(handles.Plot_IR,timeIR,selIR);
    end
    %Si ya hay resultados:
    if (ResultsObtained)
        ShowGridResults(hObject,handles,selpos,source)
    end

function List_Positions_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function Edit_XPix_Callback(hObject, eventdata, handles)

function Edit_XPix_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function Edit_YPix_Callback(hObject, eventdata, handles)

function Edit_YPix_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function Edit_XMeters_Callback(hObject, eventdata, handles)

function Edit_XMeters_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function Edit_YMeters_Callback(hObject, eventdata, handles)

function Edit_YMeters_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function Btn_DeletePos_Callback(hObject, eventdata, handles)
    global Positions npositions
    selpos = get(handles.List_Positions,'Value');
    list = get(handles.List_Positions,'String');
    list = [list(1:selpos-1);list(selpos+1:end)];
    Positions = [Positions(1:selpos-1,:,:);Positions(selpos+1:end,:,:)];
    npositions = npositions-1;
    if(get(handles.List_Positions,'Value')==length(list)+1)
        set(handles.List_Positions,'Value',selpos-1)
    end
    set(handles.List_Positions,'String',list);
    PlotPositions(handles,hObject,0);

function Btn_PickPos_Callback(hObject, eventdata, handles)

    global picked npositions
    %Habilita la seleccion por click en la imagen del plano e incrementa en 1
    %el numero de posiciones.
    if picked == 0
        picked = 1;
        set(handles.imageHandle,'ButtonDownFcn',{@MapClickCallback,handles,hObject});
        guidata(hObject,handles);
        npositions = npositions + 1;
    end

function Edit_PositionName_Callback(hObject, eventdata, handles)

function Edit_PositionName_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function Btn_Calculate_Callback(hObject, eventdata, handles)

    global npositions Positions ResultsObtained NSources
    calcwait = waitbar(0,'Preparing Smoke Bomb.');
    
    %Obtención de parámetros de entrada.
    
    isoctave = get(handles.IsOct,'Value');
    isthirdoctave = get(handles.Is3Oct,'Value');
    NCIsNone = get(handles.NCIsNone,'Value');
    NCIsLundeby = get(handles.NCIsLundeby,'Value');
    NCIsPepino = get(handles.NCIsPepino,'Value');
    
    if (isoctave)
        octave = 1;
    elseif (isthirdoctave)
        octave = 3;
    end

    if (NCIsNone)
        NC = 0;
    elseif (NCIsLundeby)
        NC = 1;
    elseif (NCIsPepino)
        NC = 2;
    end

    %Settings:
    threshold = -20;
    %Calculo de parámetros para cada posición:
    AllResults = cell(NSources + 1,npositions); %2 para las fuentes y 1 para la media de ambas
    fc = 0;
    for i=1:npositions
        waitbar(i/npositions)
        Type = cell2mat(Positions(i,5,1));
        if (Type == 1) %Esto es monoaural:
            for k=1:NSources          
                %Las fs de todos los audios deben ser iguales!
                [IR,fs] = audioread(char(Positions(i,6,k)));
                CalRMS = cell2mat(Positions(i,12,k));  
                [WN,~] = audioread(char(Positions(i,7,k)));
                [AN1,~] = audioread(char(Positions(i,8,k)));
                [AN2,~] = audioread(char(Positions(i,9,k)));
                [AN3,~] = audioread(char(Positions(i,10,k)));
                [AN4,~] = audioread(char(Positions(i,11,k)));        
                disp(Positions(i,6,k));
                [Parameters,EchoSpeech,EchoMusic,fc,FS] = CalculateParametersMono(IR',WN,fs,octave,NC,threshold,CalRMS); %Calcula todos los parametros acusticos y los mete en un cell array  
                [TauMins,DeltaTau] = CalculateAnechoics(AN1,AN2,AN3,AN4,handles.TauMinsAnechoics,fs);                 
                Taus = zeros(4,length(fc)+4);
                DeltaTaus = zeros(4,length(fc)+4);
                Taus(1,1) = TauMins(1);
                Taus(2,1) = TauMins(2);
                Taus(3,1) = TauMins(3);
                Taus(4,1) = TauMins(4);
                DeltaTaus(1,1) = DeltaTau(1);
                DeltaTaus(2,1) = DeltaTau(2);
                DeltaTaus(3,1) = DeltaTau(3);
                DeltaTaus(4,1) = DeltaTau(4);
                Parameters = vertcat(Parameters,Taus,DeltaTaus);
                AllResults{k,i} = Parameters;                      
                Echoes{1,i,k} = EchoSpeech;
                Echoes{2,i,k} = EchoMusic;

            end
            %Computo de la media entre posiciones de fuente:
        MeanParametersi = 0;
            for l = 1:NSources            
                MeanParametersi = MeanParametersi + AllResults{l,i};
            end
            AllResults{NSources+1,i} = MeanParametersi/NSources;  
        elseif (Type == 2) %Binaural:
            
            for k = 1:NSources
                
                [IR,fs] = audioread(char(Positions(i,6,k)));
                [AN1,~] = audioread(char(Positions(i,8,k)));
                [AN2,~] = audioread(char(Positions(i,9,k)));
                [AN3,~] = audioread(char(Positions(i,10,k)));
                [AN4,~] = audioread(char(Positions(i,11,k))); 
                disp(Positions(i,6,k));
                [ParametersDH,fc,Fs] = CalculateParametersStereo(IR,AN1,AN2,AN3,AN4,fs,handles.RMSCalDH,octave,NC,-20);
                AllResults{k,i} = ParametersDH;
                
            end
            %Computo de la media entre posiciones de fuente:
        MeanParametersi = 0;
            for l = 1:NSources               
                MeanParametersi = MeanParametersi + AllResults{l,i};
            end
            AllResults{NSources+1,i} = MeanParametersi/NSources;
        elseif (Type == 3) %SF
            
            for k = 1:NSources
                [IR,fs] = audioread(char(Positions(i,6,k)));
                disp(Positions(i,6,k));
                [ParametersSF] = CalculateSF(IR,fs,octave,NC,threshold);
                AllResults{k,i} = ParametersSF;
            end
            
        end
        %Computo de la media entre posiciones de fuente:
        MeanParametersi = 0;
            for l = 1:NSources               
                MeanParametersi = MeanParametersi + AllResults{l,i};
            end
            AllResults{NSources+1,i} = MeanParametersi/NSources;
    end
    save('results.mat','AllResults')
    handles.AllResults = AllResults;
    handles.fc = fc;
    handles.FS = FS;
    handles.echoes = Echoes;

    ResultsObtained = 1;

    %Seteo de resultados:

    ParametersString = {'C50';'C80';'D50';'EDT';'T20';'T30';'STI';'ALCons %';'Echo Music Max';'Echo Music Time';'Echo Speech Max';'Echo Speech Time';'DR';'STEarly1';'STEarly2';'STLate';'RR160';'SPL';'TauE Min (1)';'TauE Min (2)';'TauE Min (3)';'TauE Min (4)';'DeltaTauE Min (1)';'DeltaTauE Min (2)';'DeltaTauE Min (3)';'DeltaTauE Min (4)';'IACCe';'IACCl';'IACCa';'IACC Ando (1)';'IACC Ando (2)';'IACC Ando (3)';'IACC Ando (4)';'LFearly';'LFlate'};
    handles.strparam = ParametersString;

    
    set(handles.Combo_Parameter,'String',ParametersString)
    colstr(1) = cellstr('Global');
    colstr(2) = cellstr('IR 10 ms');
    colstr(3) = cellstr('IR 100 ms');
    colstr(4) = cellstr('IR 350 ms');
    col2str(1) = cellstr('Global');
    
    for i = 5:4+length(fc)
        colstr(i) = cellstr(num2str(fc(i-4)));
    end
    for i = 2:1+length(fc)
        col2str(i) = cellstr(num2str(fc(i-1)));
    end    
    handles.colstr = colstr;
    handles.col2str = col2str;
    set(handles.Combo_Band,'String',colstr)
    set(handles.Combo_Parameter,'Value',1)
    set(handles.Combo_Band,'Value',1)
    close(calcwait)
    PlotMap(handles,hObject,1,1,1,1)
    guidata(hObject,handles);
    ShowGridResults(hObject,handles,1,1)

function ShowGridResults(hObject,handles,selectedpos,selectedsource)
    global Positions

    AllResults = handles.AllResults;
    colstr = handles.colstr;
    col2str = handles.col2str;
    
    RowStringMON = {'C50';'C80';'D50';'EDT';'T20';'T30';'STI';'ALCons %';'Echo Music Max';'Echo Music Time';'Echo Speech Max';'Echo Speech Time';'D/R';'STEarly1';'STEarly2';'STLate';'RR160';'SPL';'TauE Min (1)';'TauE Min (2)';'TauE Min (3)';'TauE Min (4)';'DeltaTauE Min (1)';'DeltaTauE Min (2)';'DeltaTauE Min (3)';'DeltaTauE Min (4)'};
    RowStringBIN = {'IACCe';'IACCl';'IACCa';'IACC Ando (1)';'IACC Ando (2)';'IACC Ando (3)';'IACC Ando (4)'};
    RowStringSF = {'LFearly','LFlate'};
    
    Type = cell2mat(Positions(selectedpos,5,selectedsource));
    switch Type
        case 1
            set(handles.Grid_Results,'RowName',RowStringMON)
            set(handles.Grid_Results,'ColumnName',colstr);
        case 2
            set(handles.Grid_Results,'RowName',RowStringBIN)
            set(handles.Grid_Results,'ColumnName',col2str);
        case 3
            set(handles.Grid_Results,'RowName',RowStringSF)
            set(handles.Grid_Results,'ColumnName',colstr);
    end
    
    set(handles.Grid_Results,'Data',AllResults{selectedsource,selectedpos})
    

function mnu_File_Callback(hObject, eventdata, handles)

function mnu_Tools_Callback(hObject, eventdata, handles)

function mnu_About_Callback(hObject, eventdata, handles)
    msgbox(sprintf(['Acuntref v1.0\n' '\n' 'Software developed by Leonardo Pepino\n' '\n' 'Copyright 2018\n' '\n' 'If you feel this application helped you in your final work, please consider donating your final note multiplied by the number of measured positions or an equivalent beer quantity.']),'About')

function mnu_Convolver_Callback(hObject, eventdata, handles)
GUI_BatchConvolver();

function mnu_Cutter_Callback(hObject, eventdata, handles)
GUI_BatchCutter();

function mnu_Open_Callback(hObject, eventdata, handles)

function mnu_SavePos_Callback(hObject, eventdata, handles)

    global Positions

    floorplan = handles.floorplan;
    polyregion = handles.polyregion;
    xvert = handles.xvert;
    yvert = handles.yvert;
    imagesize = handles.imagesize;
    alpha = handles.alpha;
    [FileName,PathName] = uiputfile('*.mat','Save Positions As');
    save(strcat(PathName,'\',FileName),'Positions','floorplan','polyregion','xvert','yvert','imagesize','alpha');
    if isfield(handles,'RMSCalDH')
        RMSCalDH = handles.RMSCalDH;
        save(strcat(PathName,'\',FileName),'RMSCalDH','-append');
    end
    if isfield(handles,'TauMinsAnechoics')
        TauMinsAnechoics = handles.TauMinsAnechoics;   
        save(strcat(PathName,'\',FileName),'TauMinsAnechoics','-append');
    end  
    if isfield(handles,'AllResults')
        AllResults = handles.AllResults;
        fc = handles.fc;
        FS = handles.FS;
        %echoes = handles.echoes;
        colstr = handles.colstr;
        col2str = handles.col2str;
        strparam = handles.strparam;  
        save(strcat(PathName,'\',FileName),'AllResults','fc','FS','colstr','col2str','strparam','-append');
    end    
    

function mnu_SavePrj_Callback(hObject, eventdata, handles)

function mnu_ExportMap_Callback(hObject, eventdata, handles)
    uiwait(msgbox('This will export alll combinations of map, thats a lot of data so if you proceed you will need space in your disk and time. Everything will be saved in an ordered way in the folder you select','modal'))    
    folder_name = uigetdir();
    nparams = length(handles.strparam);
    nsources = 3;
    mapwait = waitbar(0,'Smoke some weed.');
    for i=1:nparams
        waitbar(i/nparams);
        mkdir(folder_name,char(handles.strparam(i))) %Crea la carpeta del parametro:
        for j=1:nsources
            mkdir(strcat(folder_name,'\',char(handles.strparam(i))),strcat('Source ',num2str(j))) %Crea las carpeta de cada fuente:
            
            if (i==7||i==8||(i>18&&i<27)||i>29&&i<34) %casos globales
                nbands = 4;
            else
                if(i<27||i>33)
                    nbands = length(handles.colstr);
                    col2use = handles.colstr;
                else
                    nbands = length(handles.col2str);
                    col2use = handles.col2str;
                end
            end
            
            for k=1:nbands
                PlotMap(handles,hObject,k,i,j,0);
                movefile('FloorMap.jpg',strcat(folder_name,'\',char(handles.strparam(i)),'\','Source ',num2str(j),'\',char(col2use(k)),'.jpg'));
            end        
        end        
    end
    close(mapwait)
    

function mnu_ExportXLS_Callback(hObject, eventdata, handles)
    
    global npositions Positions
    nparams = length(handles.strparam);
    AllResults = handles.AllResults;
    xlswait = waitbar(0,'Wait some seconds please');
    for source = 1:3
        
        waitbar(source/3)
        
        for param=1:nparams
            i = 1;
            %Cols = 0;
            XLS = 0;
            
            if (param==7||param==8||(param>18&&param<27)||param>29&&param<34) %casos globales
                nbands = 4;
            else
                if(param<27||param>33)
                    nbands = length(handles.colstr);
                    col2use = handles.colstr;
                else
                    nbands = length(handles.col2str);
                    col2use = handles.col2str;
                end
            end
            clear 'Cols'
            for band = 1:nbands
                Cols{1,band} = char(col2use(band));
            end
            
            for pos = 1:npositions
                
                Type = cell2mat(Positions(pos,5,1));
                if (param>0 && param<27 && Type == 1)   
                    
                    for band = 1:nbands
                        XLS(i,band) = AllResults{source,pos}(param,band);
                    end
                    posstr{i} = num2str(pos);
                    i = i+1;
                elseif (param>26 && param<34 && Type == 2)
                    
                    for band = 1:nbands
                        
                        XLS(i,band) = AllResults{source,pos}(param-26,band);
                    end
                    posstr{i} = num2str(pos);
                    i = i+1;
                    
                elseif (param>33 && Type == 3)
                    
                    for band = 1:nbands
                        XLS(i,band) = AllResults{source,pos}(param-33,band);
                    end
                    posstr{i} = num2str(pos);
                    i = i+1;
                end
            end
            posstr = posstr(1:i-1);
            xlswrite(strcat('Source ',num2str(source)),Cols,char(handles.strparam(param)),'B1');
            xlswrite(strcat('Source ',num2str(source)),XLS,char(handles.strparam(param)),'B2');
            xlswrite(strcat('Source ',num2str(source)),posstr',char(handles.strparam(param)),'A2');
        end      
    end
    close(xlswait);
    
function mnu_FilePos_Callback(hObject, eventdata, handles)

    global Positions npositions ResultsObtained
    
    [filename, pathname] = uigetfile('*.mat','Select positions file');
    S = load(strcat(pathname,'\',filename));
    Positions = S.Positions;
    handles.alpha = S.alpha;
    handles.floorplan = S.floorplan;
    handles.imagesize = S.imagesize;
    handles.polyregion = S.polyregion;
    handles.xvert = S.xvert;
    handles.yvert = S.yvert;
    
    if isfield(S,'RMSCalDH')
        handles.RMSCalDH = S.RMSCalDH;
    end
    
    if isfield(S,'TauMinsAnechoics')
        handles.TauMinsAnechoics = S.TauMinsAnechoics;
    end
    
    if isfield(S,'AllResults')
        handles.AllResults = S.AllResults;
        handles.fc = S.fc;
        handles.FS = S.FS;
        %handles.echoes = S.echoes;
        handles.colstr = S.colstr;
        handles.col2str = S.col2str;
        handles.strparam = S.strparam;    
        set(handles.Combo_Parameter,'String',handles.strparam)
        set(handles.Combo_Band,'String',handles.colstr)
        set(handles.Combo_Parameter,'Value',1)
        set(handles.Combo_Band,'Value',1)
        ResultsObtained = 1;
    end
    
    axes(handles.Plot_Map);
    handles.imageHandle = imshow(S.floorplan);
    npos = size(Positions);
    npositions = npos(1);
    for i=1:npositions
        listpositions(i) = Positions(i,1,1);
    end
    set(handles.List_Positions,'String',listpositions);

    
    guidata(hObject,handles);
    PlotPositions(handles,hObject,0);
    
    if isfield(S,'AllResults')
        ShowGridResults(hObject,handles,1,1)
    end
    

function mnu_OpenPrj_Callback(hObject, eventdata, handles)

function Btn_LoadPlan_Callback(hObject, eventdata, handles)

    global picked
    %Abre la imagen y la muestra:
    [filename, pathname] = uigetfile({'*.jpg';'*.png';'*.tiff';'*.*'},'Select an Image');
    set(gcf,'CurrentAxes',handles.Plot_Map)
    floorplan = imread(strcat(pathname,filename));
    handles.floorplan = floorplan;
    handles.imageHandle = imshow(floorplan);

    %Definir area de mapeo:

    uiwait(msgbox('Nice smoke. Now draw a polygon defining area to map. When closed, double click inside it to finish.','Define Mapping Area','modal'));
    [maparea,xvert,yvert] = roipoly(floorplan);
    handles.polyregion = maparea;
    handles.xvert = xvert;
    handles.yvert = yvert;
    handles.imagesize = size(maparea);

    alpha = ones(handles.imagesize(1),handles.imagesize(2));
     for x = 1:handles.imagesize(1)
         for y = 1:handles.imagesize(2)
             if (maparea(x,y))
                 alpha(x,y) = 0.2;
             end
         end
     end

    handles.imageHandle = imshow(floorplan);
    handles.alpha = alpha;
    % set(handles.imageHandle,'AlphaData',alpha);
    uiwait(msgbox('Well Done. Now create listener points and assign IRs. This can be a paja so once finished consider saving positions from File Menu. To create listener points, click on Pick Position and then click over the floor plan in the desired position. Repeat the process of clicking pick position and floor plan. Once finished, go to File menu and load monoaural calibration files.','Mapping Area Defined','modal'))

    guidata(hObject,handles);
    picked = 0;

function MapClickCallback(objectHandle,~,handles,hObject)

    %Evento para dibujar puntos al hacer click
    global picked
    if (picked==1)
        axesHandle = get(objectHandle,'Parent');
        coordinates = get(axesHandle,'CurrentPoint');
        pixcoordinates = coordinates(1,1:2);
        hold(axesHandle,'on')
        plot(axesHandle,pixcoordinates(1,1),pixcoordinates(1,2),'o','MarkerSize',8,'MarkerFaceColor','k');
        picked = 0;    
        AddNewPosition(handles,hObject,pixcoordinates);
        guidata(hObject,handles);
    end
       
function AddNewPosition(handles,hObject,pixcoordinates)

    global npositions Positions counter NSources
    counter = counter+1;
    for i = 1:NSources
        Positions{npositions,1,i} = num2str(counter);
        Positions{npositions,2,i} = pixcoordinates(1,1);
        Positions{npositions,3,i} = pixcoordinates(1,2);
    end
    set(handles.List_Positions,'String',Positions(:,1,1));
    set(handles.List_Positions,'Value',npositions);
    guidata(hObject,handles);
    
function PlotPositions(handles,hObject,highlightedi)

    global Positions
    npos = size(Positions);
    cla(handles.Plot_Map);
    axes(handles.Plot_Map)
    handles.imageHandle = imshow(handles.floorplan);
    hold(handles.Plot_Map,'on')
    for i=1:npos(1)
        plot(handles.Plot_Map,cell2mat(Positions(i,2,1)),cell2mat(Positions(i,3,1)),'o','MarkerSize',8,'MarkerFaceColor','k');
        hold(handles.Plot_Map,'on')
        
        if i==highlightedi
            
            plot(handles.Plot_Map,cell2mat(Positions(i,2,1)),cell2mat(Positions(i,3,1)),'o','MarkerSize',8,'MarkerFaceColor','g');
            hold(handles.Plot_Map,'on')
            
        end
        
    end
    guidata(hObject,handles);

function PlotMap(handles,hObject,band,param,source,show)
    % 1 - 26 -> MON
    % 27 - 33 -> BIN
    % 34 -> SF
    global Positions npositions
    if show
        mapwait = waitbar(0,'Generating smoke final attack. This could take some seconds maybe minutes. 500 Hz. Coff coff');
    end
    AllResults = handles.AllResults;
    %band = get(handles.Combo_Band,'Value');
    %param = get(handles.Combo_Parameter,'Value');
    %source = get(handles.Combo_SourceMap,'Value');
    k = 1;
    
    if param < 27
        for i=1:npositions
            Typei = cell2mat(Positions(i,5,1));
            if  Typei == 1
                Parametertomap(k) = AllResults{source,i}(param,band);
                Xtomap(k) = cell2mat(Positions(i,2,1));
                Ytomap(k) = cell2mat(Positions(i,3,1));
                k = k+1;
            end
        end
    elseif (param > 26 && param <34)
        for i=1:npositions
            Typei = cell2mat(Positions(i,5,1));
            if  Typei == 2
                Parametertomap(k) = AllResults{source,i}(param-26,band);
                Xtomap(k) = cell2mat(Positions(i,2,1));
                Ytomap(k) = cell2mat(Positions(i,3,1));
                k = k+1;
            end
        end
    elseif param == 34
        for i=1:npositions
            Typei = cell2mat(Positions(i,5,1));
            if  Typei == 3
                Parametertomap(k) = AllResults{source,i}(param-33,band);
                Xtomap(k) = cell2mat(Positions(i,2,1));
                Ytomap(k) = cell2mat(Positions(i,3,1));
                k = k+1;
            end
        end
    end
    if show
        waitbar(1/4)
    end
    [xf,yf,map] = mapping(Xtomap',Ytomap',Parametertomap,handles.xvert,handles.yvert,handles.imagesize);
    if show
        waitbar(2/4)
    end
    %Se crea una imagen del mapeo con el mismo tamaño de la imagen del plano y
    %se guarda en archivo temporal:
    mapfig = figure('rend','opengl','Visible','off');
    xf = 1:handles.imagesize(2);
    yf = 1:handles.imagesize(1);
    [hmap,cobj] = contourf(xf,yf,map',200);
    cobj.LineStyle = 'none';
    colorstr = get(handles.Combo_Colormap,'String');
    selectedcolor = get(handles.Combo_Colormap,'Value');
    contourcmap(char(colorstr(selectedcolor)));
    caxis manual;
    caxis([min(Parametertomap),max(Parametertomap)])
    cbar = colorbar;
    cbar.Units = 'pixels';
    cbarposition = get(cbar,'Position');
    cbarwidth = round(cbarposition(3));
    cbar.Position = [handles.imagesize(2)+cbarwidth/2,10,cbarwidth,handles.imagesize(1)-20];
    mapfig.Position = [0 0 handles.imagesize(2) + cbarwidth*4 handles.imagesize(1)];
    ax = gca;
    dpi = 100;
    ax.XTick = [];
    ax.YTick = [];
    ax.YMinorGrid = 'off';
    ax.YMinorTick = 'off';
    ax.XMinorGrid = 'off';
    ax.XMinorTick = 'off';
    ax.Units = 'pixels';
    ax.OuterPosition = [0 0 handles.imagesize(2) handles.imagesize(1)];
    ax.Position = [0 0 handles.imagesize(2) handles.imagesize(1)];
    set(mapfig,'PaperUnits','inches','PaperPosition',[0 0 ((handles.imagesize(2)) + cbarwidth*4)/dpi handles.imagesize(1)/dpi])
    if show
        waitbar(3/4)
    end
    print(mapfig,'mapeoo3.jpg','-djpeg','-r100')

    wholefig = figure('rend','opengl','Visible','off');
    wholefig.Position = [0 0 handles.imagesize(2) + cbarwidth*4 handles.imagesize(1)];
    mapimage = imread('mapeoo3.jpg');
    hmap = imshow(mapimage);
    hold on
    him = imshow(handles.floorplan);
    alpha = handles.alpha;
    set(him,'AlphaData',alpha); %Setea transparencia en area de mapeo.
    set(wholefig,'PaperUnits','inches','PaperPosition',[0 0 ((handles.imagesize(2)) + cbarwidth*4)/dpi handles.imagesize(1)/dpi])
    figtitle = strcat(handles.strparam(param),' - ',handles.colstr(band));
    
    if (band ~= 1)
        figtitle = strcat(figtitle,' Hz');
    end
    
    if (source ~= 3)
        figtitle = strcat(figtitle,' Source Position: ',num2str(source));
    else
        figtitle = strcat(figtitle,' Mean of Source Positions');
    end
    
    title(figtitle)
    print(wholefig,'FloorMap.jpg','-djpeg','-r150')

    if (show)
        floormap = imread('FloorMap.jpg');
        cla(handles.Plot_Map);
        axes(handles.Plot_Map);
        handles.imageHandle = imshow(floormap);
    
        waitbar(4/4);
        close(mapwait);
    end
    
function Btn_IRLoad_Callback(hObject, eventdata, handles)

    [filename, pathname] = uigetfile({'*.wav';'*.mp3'},'Select an impulse response audio file');

    wholepath = strcat(pathname,filename);
    [selIR,fs] = audioread(wholepath);

    set(handles.Edit_IRFile,'String',wholepath);

    timeIR = (0:length(selIR)-1)/fs;
    plot(handles.Plot_IR,timeIR,selIR);
    axis(handles.Plot_IR,[0 timeIR(length(timeIR)) -1 1]);

function Edit_IRFile_Callback(hObject, eventdata, handles)

function Edit_IRFile_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function Btn_WNLoad_Callback(hObject, eventdata, handles)

    [filename, pathname] = uigetfile({'*.wav';'*.mp3'},'Select a white noise audio file');
    wholepath = strcat(pathname,filename);

    set(handles.Edit_WNFile,'String',wholepath);

function Edit_WNFile_Callback(hObject, eventdata, handles)

function Edit_WNFile_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function Btn_AN1Load_Callback(hObject, eventdata, handles)

    [filename, pathname] = uigetfile({'*.wav';'*.mp3'},'Select Anechoic 1 audio file');
    wholepath = strcat(pathname,filename);

    set(handles.Edit_AN1File,'String',wholepath);

function Edit_AN1File_Callback(hObject, eventdata, handles)

function Edit_AN1File_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function Btn_AN2Load_Callback(hObject, eventdata, handles)

    [filename, pathname] = uigetfile({'*.wav';'*.mp3'},'Select Anechoic 2 audio file');
    wholepath = strcat(pathname,filename);

    set(handles.Edit_AN2File,'String',wholepath);

function Edit_AN2File_Callback(hObject, eventdata, handles)

function Edit_AN2File_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function Btn_AN3Load_Callback(hObject, eventdata, handles)

    [filename, pathname] = uigetfile({'*.wav';'*.mp3'},'Select Anechoic 3 audio file');
    wholepath = strcat(pathname,filename);

    set(handles.Edit_AN3File,'String',wholepath);

function Edit_AN3File_Callback(hObject, eventdata, handles)

function Edit_AN3File_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function Btn_AN4Load_Callback(hObject, eventdata, handles)

    [filename, pathname] = uigetfile({'*.wav';'*.mp3'},'Select Anechoic 4 audio file');
    wholepath = strcat(pathname,filename);

    set(handles.Edit_AN4File,'String',wholepath);

function Edit_AN4File_Callback(hObject, eventdata, handles)

function Edit_AN4File_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function Combo_PosType_Callback(hObject, eventdata, handles)

function Combo_PosType_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function Edit_MicNumber_Callback(hObject, eventdata, handles)

function Edit_MicNumber_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function Btn_Apply_Callback(hObject, eventdata, handles)

    global Positions

    selpos = get(handles.List_Positions,'Value');
    source = get(handles.Combo_SourceSet,'Value'); %Despues ponerlo en funcion de un elemento de la gui
    Positions{selpos,4,source} = str2num(get(handles.Edit_MicNumber,'String'));
    Positions{selpos,5,source} = get(handles.Combo_PosType,'Value');
    Positions{selpos,6,source} = get(handles.Edit_IRFile,'String');
    Positions{selpos,7,source} = get(handles.Edit_WNFile,'String');
    Positions{selpos,8,source} = get(handles.Edit_AN1File,'String');
    Positions{selpos,9,source} = get(handles.Edit_AN2File,'String');
    Positions{selpos,10,source} = get(handles.Edit_AN3File,'String');
    Positions{selpos,11,source} = get(handles.Edit_AN4File,'String');

function Btn_CalcPix_Callback(hObject, eventdata, handles)

    if (handles.ScaleDefined == 0)
        [x, y] = getpts(handles.Plot_Map);
        pixdistance = sqrt((x(2)-x(1)).^2+(y(2)-y(1)).^2);    
        prompt = {'Enter measured distance:'};
        dlg_title = 'Floor Plan Scale Setting';
        num_lines = 1;
        defaultans = {'20'};
        mdistance = str2num(char((inputdlg(prompt,dlg_title,num_lines,defaultans))));
        handles.scale = mdistance/pixdistance;
        handles.ScaleDefined = 1;
        Xpix = str2num(get(handles.Edit_XPix,'String'));
        Ypix = str2num(get(handles.Edit_YPix,'String'));
        XMeters = Xpix*handles.scale;
        YMeters =(handles.imagesize(1)-Ypix)*handles.scale;
        set(handles.Edit_XMeters,'String',num2str(XMeters));
        set(handles.Edit_YMeters,'String',num2str(YMeters));   
    end
    guidata(hObject,handles);

function Grid_Results_CellSelectionCallback(hObject, eventdata, handles)
    source = get(handles.Combo_SourceMap,'Value');
    selection = eventdata.Indices;
    selrow = selection(1);
    selcol = selection(2);
    pos = get(handles.List_Positions,'Value');
    AllResults = handles.AllResults;
    selparam = AllResults{source,pos}(selrow,:);
    plot(handles.Plot_Parameter,selparam);
 
    
function mnu_BatchLoad_Callback(hObject, eventdata, handles)

    global Positions
    [filename_audios, pathname_audios] = uigetfile({'*.wav','WAV-files (*.wav)'},'Pick audio files','MultiSelect', 'on');

    for i=1:length(filename_audios)
        nameswav = strsplit(char(filename_audios(i)),'.');
        nameswav = strsplit(char(nameswav(1)),'_'); %Separo el nombre en bloques para interpretar a donde va:
        positionstr = char(nameswav(4));
        positionidx = str2num(positionstr(2:end));
        source = char(nameswav(3));
        source = str2num(source(2:end));
        switch (char(nameswav(2)))
            case 'AN1'
                Positions{positionidx,8,source} = strcat(pathname_audios,'\',char(filename_audios(1,i)));
            case 'AN2'
                Positions{positionidx,9,source} = strcat(pathname_audios,'\',char(filename_audios(1,i)));
            case 'AN3'
                Positions{positionidx,10,source} = strcat(pathname_audios,'\',char(filename_audios(1,i)));
            case 'AN4'
                Positions{positionidx,11,source} = strcat(pathname_audios,'\',char(filename_audios(1,i)));
            case 'SWEEP'
                Positions{positionidx,6,source} = strcat(pathname_audios,'\',char(filename_audios(1,i)));
            case 'WHITE'
                Positions{positionidx,7,source} = strcat(pathname_audios,'\',char(filename_audios(1,i)));
        end
        micstr = char(nameswav(1));
        mictype = micstr(1:2);
        switch mictype
            case 'EW'
                Positions{positionidx,5,source} = 1;
                micnumber = str2num(char(micstr(3:end)));
                Positions{positionidx,4,source} = micnumber;
                RMSCal = handles.RMSCal(micnumber);
                Positions{positionidx,12,source} = RMSCal;
            case 'DH'
                Positions{positionidx,5,source} = 2;
                Positions{positionidx,4,source} = -1;
                RMSCal = handles.RMSCalDH;
                Positions{positionidx,12,source} = RMSCal;
            case 'SF'
                Positions{positionidx,5,source} = 3;
                Positions{positionidx,4,source} = -1;
                Positions{positionidx,12,source} = -1;
        end       
    end
    
    msgbox('Well done, now you are ready for the final fight! If you want to modify any file association with an specific position, just select the position in the list and click on the corresponding button on File Panel. If you feel you are ready just select a Noise correction method... coff coff Pepinos one not Lundeby... and the third octave filtering as Smoke likes most. Finally click on the magic button Calculate and wait not too much time I promise.','You are almost there')
    guidata(hObject,handles);
    
function mnu_LoadCal_Callback(hObject, eventdata, handles)

    %Acá abrir calibraciones, con algun codigo de nombre, calcular RMSs de
    %referencia y meterlos en handles.RMSCal
    [filename_audios, pathname_audios] = uigetfile({'*.wav','WAV-files (*.wav)'},'Pick calibration files','MultiSelect', 'on');
    for i=1:length(filename_audios)
        micname = strsplit(char(filename_audios(i)),'.');
        micnumber = str2num(char(micname(1)));
        [xcal,fs] = audioread(strcat(pathname_audios,'\',char(filename_audios(1,i))));
        calRMS = CalcRMS(xcal,fs,1);
        RMSCal(micnumber) = calRMS;
    end
    handles.RMSCal = RMSCal;
    
    msgbox('Well done, all Earthwork microphones are calibrated with laboratory precision. Now load those smoky and irritating anechoic files you played during the smoke session over and over again. Note that Sato code is dangerous and will throw lots of warnings. Later I will deactivate them, now just ignore them.','You are a winner')
    guidata(hObject,handles);

function mnu_LoadANs_Callback(hObject, eventdata, handles)

    %Aca abrir anecoicos, con algun código de nombre, calcular tauemins y
    %meterlos en handles.TauMinsAnechoics
    [~, pathname_audios] = uigetfile({'*.wav','WAV-files (*.wav)'},'Pick anechoic files','MultiSelect', 'on');
    [AN1,fs] = audioread(strcat(pathname_audios,'\','AN1.wav'));
    [AN2,~] = audioread(strcat(pathname_audios,'\','AN2.wav'));
    [AN3,~] = audioread(strcat(pathname_audios,'\','AN3.wav'));
    [AN4,~] = audioread(strcat(pathname_audios,'\','AN4.wav'));
    Tau(1) = TauEMin(AN1,fs);
    Tau(2) = TauEMin(AN2,fs);
    Tau(3) = TauEMin(AN3,fs);
    Tau(4) = TauEMin(AN4,fs);
    handles.TauMinsAnechoics = Tau;
    
    msgbox('You are almost there. Now calibrate binaural files. You have to first open Left channel calibration file, and then right one','Developing love to slave work')
     
    guidata(hObject,handles);

function mnu_CalibrateBinaural_Callback(hObject, eventdata, handles)
    
    uiwait(msgbox('Select left channel calibration file','Binaural Calibration','modal'));
    [fileL, pathL] = uigetfile({'*.wav','WAV-files (*.wav)'},'Pick Left Channel Calibration File');
    uiwait(msgbox('Now select right channel calibration file','Binaural Calibration','modal'));
    [fileR, pathR] = uigetfile({'*.wav','WAV-files (*.wav)'},'Pick Right Channel Calibration File');
    [CalL,fsL] = audioread(strcat(pathL,'\',fileL));
    [CalR,fsR] = audioread(strcat(pathR,'\',fileR));
    
    RMSL = CalcRMS(CalL,fsL,1);
    RMSR = CalcRMS(CalR,fsR,1);
    handles.RMSCalDH = RMSR/RMSL;
    
    msgbox('Well done. The japanese duet, Ando and the Satos is proud of you. You are one step closer to the final boss. Prepare for the fight and go to Batch Load and select aaaaallllll of the measured audios. If you dont want me to go angry the filenames should have the next format : (microphone)_(content)_(source)_(position).wav. Microphones can be EW(mic number) (Earthworks),SF(Soundfield),DH(Ale Head), and should be followed by the mic number. Content should be AN(number), SWEEP or WHITE. Finally, source is F(number), and position P(number). For example: EW7_AN3_F1_P7.wav corresponds to the anechoic excerpt 3 recorded with Earthwork mic number 7 in position 7 with source in position 1. Understood? Finally note that dummy head files should be stereo, and SF also, with 8 figure in left channel and omni in right.','Please read this carefully')
    guidata(hObject,handles);
    
function Combo_SourceMap_Callback(hObject, eventdata, handles)    
    
function Combo_SourceMap_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Combo_SourceSet_Callback(hObject, eventdata, handles)

    List_Positions_Callback(hObject, eventdata, handles)
       
function Combo_SourceSet_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Btn_DrawMap.
function Btn_DrawMap_Callback(hObject, eventdata, handles)
    
    band = get(handles.Combo_Band,'Value');
    param = get(handles.Combo_Parameter,'Value');
    source = get(handles.Combo_SourceMap,'Value');
    
    PlotMap(handles,hObject,band,param,source,1)

% --- Does nothing, but there are not more warnings now!
function mnu_Export_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function mnu_ImportXLS_Callback(hObject, eventdata, handles)
    
    uiwait(msgbox('Select modified xls files','Import XLS','modal'));
    
    [filename_XLS, pathname_XLS] = uigetfile({'*.xls;*.xlsx'},'Microsoft Excel Spreadsheets','MultiSelect', 'on');
    nparams = 35;
    xlsiwait = waitbar(0,'Wait some seconds please');
    for source = 1:length(filename_XLS)
        waitbar(source/3);
        sourcenumber = strsplit(char(filename_XLS(source)),'e');
        sourcenumber = strsplit(char(sourcenumber(2)),'.');
        sourcenumber = str2num(char(sourcenumber(1)));
        filepath = strcat(pathname_XLS,'\',char(filename_XLS(1,source)));
        
        for param = 1:nparams
            
            k = param+1;
            
            Parami = xlsread(filepath,k);
            sizexls = size(Parami);
            if (param==7||param==8||(param>18&&param<27)) %casos globales
                nbands = 4;
                a = 1;
            else
                if (param>29&&param<34) %Solo global
                    a = 2;
                    nbands = 1;
                else
                    nbands = sizexls(2)-1;
                    a = 2;
                end
            end 
            
            if (param>26 && param<34)
                param = param-26;
            elseif (param>33)
                param = param - 33;
            end
            
            for pos = a:sizexls(1)

                position = Parami(pos,1);
                for band = 2:nbands+1
                    valuecell = Parami(pos,band);
                    AllResults{source,position}(param,band-1) = valuecell;
                end
            end
        end
    end
    handles.AllResults = AllResults;
    close(xlsiwait);
    
    guidata(hObject,handles);
    ShowGridResults(hObject,handles,1,1)
