# Acountref-v-1.1
Acountref is a MATLAB software designed for the processing of acoustic measurements in a room. Room acoustic parameters such as T30, C50, and many more (including Ando parameters) are obtained for every measured position and a mapping is created. This software was developed by Leonardo Pepino in order to accelerate the processing of acoustic measurements as it takes a long time in commercial softwares to analyze many audios and organisate the data. Acoustics measurements of this type are performed in the subject Instrumentos y Mediciones Acústicas of the career Ingenieria en Sonido at UNTref university (Argentina). 

Note: Reccomended screen resolution for running this application is 1600 x 900. If Graphical interface components are not shown clearly, please test with different screen resolutions.

Instructions:

Executing GUI_Main.m opens a graphical user interface. Follow these steps:

0) You can access to a set of tools (in the Tools Menu) to improve the efficiency of your work:
    Batch Cutter: allows you to open as many audio files as you want and cut them all at the same specified positions. So one of the audio waveforms will be shown, click on Add Segment, give it a name and click in the waveform to specify the beginning and ending of the segment. Repeat as many times as you need. Then click on Segment and Save and select the output folder. All files will be cutted in the specified segments and the name of the segment will be used as a preffix in the filenames.
    
    Batch Convolver: allows to convolve many files at once. Just select the audio files, the inverse filter, the output directory and click on Go! Note: select as output directory a different folder than the audio files one to avoid overwriting and losing original data.
    
    Thanks to FedeBosio for adding user exceptions.
    
    Now the main software:

1) Press "Load Floor Plan" button to open a floor plan image of the room you measured. Over this image the mapping of parameters will be performed.

2) The image will be shown in the interface and then clicking on it will add a vertex of a polygon defining the area of the image where mapping is desired. Click and add as many vertexs as you wish and click on the first vertex to close polygon. Once closed you can drag it or just click two times inside it to finish this step.

3) Now you can add the positions you measured. Press "Pick Position" button and then click over the floor plan in the desired ubication.
  Repeat this step to add as many positions as you want.
  
4) Once you finished adding positions, please save the progress pressing Save Positions at File Menu. Now you have 2 options:
    
    a) Select each created position from the list and set manually the audio files for each position as well as the microphone number and type. (Recommended for few measured positions)
    b) Name your measurements as follows:
    
    (Mic Type and Number)_(Audio content)_(Source position)_(Measurement position).wav
    
    Mic Type can be EW (Earthworks), DH(Dummy Head), SF(Soundfield).
    Audio content can be SWEEP (Impulse response), AN[number] (Anechoic music (up to 4 different excerpts)), WHITE (White noise)
    
    For example:
    
    EW14_AN3_F1_P15.wav corresponds to the measurement taken at position 15 with Earthworks microphone number 14 during the playback of anechoic file number 3 with source position 1. 
    Dummy Head and Soundfields mics dont have number.
    
5) If you went for option b, load monoaural calibration files (which should be named as micnumber.wav -> For example: 5.wav is the calibration file of the mic 5)

6) Then load anechoic files which should be named AN[number].wav. These are the audios you played during measurements.

7) Load Dummy head calibration files (first left channel then right)

8) Load all measurement audios named as stated in 4b)

9) Now you can click on a position in the list and check if the files were assigned.

10) Select noise correction method and filtering type (Lundeby not working properly)

11) Press calculate and wait some time

12) A mapping of C50 will be displayed and results will be shown in the table. Select a position in the list to see the parameters measured there in the table. Set the parameter and desired band to show in the map and click on draw map. Each time an individual mapping is performed, the corresponding .jpeg file is saved in the root directory of the project if you need just a single mapping.

13) You can export results to an excel file from file menu. Click on it and a file for each source position will be created. It will have a sheet for each parameter and all positions results.

14) You can export all map combinations images. This will take a very long time. Every file will be ordered in folders and subfolders for a good organisation of results.

15) You can also reimport an xls file you exported but modified (maybe some values are outliers and you wish not to show them in the mapping)

16) Read the about.

Added functionalities to version 1.0:

- Now project can be saved at any progress (when you finish placing positions, when you have assigned files or when you have obtained results).
- Subjective integration windows added: 10, 100 and 350 ms. (This calculation takes time so you will have to wait significantly more to obtain results).
- XLS export bugs fixed.
- Windowing method changed in energy ratio calculations (the signal is windowed before filtering)
- XLS import option.


