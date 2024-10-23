% %% MATLAB automation of LIBS experiments
%   Theta motor controller on COM1:39400 8N1 NO FLOW CONTROL
%   !!!!! accesing the COM objects WinX32 requires elevation, Matlab must
%   be run as administrator
%   v01 First version with basic data capture
%       Real test with a modern PL (exterior) 15april2014
%   v02 creates master spectra array for further processing
%       background capture is not useful, it is better to remove ALL
%       background (emission + electronic noise + ...) by processing
%   v03 two external loops for two params
%   v05 version left by Artzai
%   v06 Camera with autofocus
%   v07 read back the spectra @from the file using an external
%   readSPE.m
%   v08 MANUAL_TRACKING: ask for a line over the video image and performs automatic shots along
%   the line in deltaR increments
%   v09 Got rid of the function, so the all the variables are available to
%   save in case of an error
%   v10 major revision so each spectrum is analyzed in real time and empty
%   spectra are discarded. problem: pulse-to-pulse delay is 0.7s ALSO
%   horizontal flipping of images. add RSD calculation
%   v11 the saturation detection is performed only for a selected
%   wavelengths, so "bad" resonance lines are not considered to discard a
%   spectrum
%   v13 includes control of laser flash lamp using an arduino (COM5) & shutter
%   control on COM6
%   v15 disconnect the firewire camera for 3 seconds using the command "2" on COM5 (arduino)
%   v16 giving a try to a "2D scanning" feature. the strategy is to move
%    almost everything to the outer loop (param2) and leave the inner loop
%    to perform the transversal scan is needed.
%   v17 frame capture moved to a function with a single location for error
%   v18 image-focus autofocus code
%   v19 pruebo la cámara china-luis
%   v20 hay unleak de memoria en el preview de la cámara, voy a poner que
%   si pasa más de 1 minuto, se pone en pausa automática
%   v21 integración del modo SOUND_STOP, será como deltaT en param1, y se
%   parará automáticamente al detectar la varianza del sonido
%   v22 hay problemas de memoria para muchos puntos, si el tamaño de
%   spectra() es > 300MB se usará un modo en el que idx2 será siempre 1 y
%   se irán guardando las iteraciones en idx2 en variables aparte
%   v23 envía twitter cuando acaba
%   v24 con AUTOGAIN en NO-AUTOCOLLECT, ajusta después de cada ráfaga,
%   v25 sin fotos en ONLY_RECORDING_PATH para que vaya más rápido
%   v26 no hace falta poner el path despues de ONLY_RECORDING_PATH
%   ToDo:  automatic calibration, X,Y,Z coordinates for every spatial point
%  atención!! para ver si se quita el error de memoria de la cámara, voy a
%  probar el comando  BCDEdit /set increaseuserva 3072
% http://es.mathworks.com/help/matlab/matlab_prog/resolving-out-of-memory-errors.html?requestedDomain=www.mathworks.com
%   v27 pruebo el método de trigger/'manual') + start(vid) para acelerar
%   los FPS de la cámara, pero lo mismo peta más, es para probar http://es.mathworks.com/help/imaq/examples/acquiring-a-single-image-in-a-loop.html
%   v28 probando el modo deltaT/2D descubrí un memoryleak. después de mucho
% probar (ver memoryLeakUltimate.m) el leak lo tiene getsnapshot() pero
% solo cuando se hace el imshow(). el caso es que hay un modo más eficiente
% de usar imshow() descrito en
% http://www.mathworks.com/matlabcentral/answers/72481-imshow-update-in-realtime-colorbar-causing-it-to-slow-down
% que voy a probar para ver si arregla lo de la memoria
%   v29  los nombres de archivos incluyen cuatro cifras 0000 porque hay
%   experimentos de más de mil puntos
%   v30 dibuja una escala graduada en la imagen si ESCALA=1
%   v31 comprueba al principio todos los puertos serie
%   v32 webcam para registrar las lecturas del medidor de energía de pulsos
%   y el bombeo si se queda en el campo de visión
%       el arduino soporta ahora comandos '3' y superior que anula esos N
%       pulsos forzando a 0.2V la línea sync
%       quito la grabación de matlabdata.mat  en cada param1 porque en
%       experimentos 2D largos alguna vez ha petado con "invalid argument"
%   v33 modo SCANLINE=2 el usuario pinta la línea
%   v34 COM2 y captura con HR2000+. A veces el ordenador va lento y se
%   pierde algún trigger del HR2000+,y se cuelga, lo he solucionado
%   capturando siempre 2 espectros menos
%   v35 medida de SNR en tiempo real e inicio de MULTI_WINDOW
%       PENDIENTE: no funciona en autocollect=1; spectros hr2000 sin hacer;
%       scaline=1 repite la línea en cada ventana; ajustar procesado rsd y
%       autogain al múltiple espectro
%       hr2000????
%   v36 de nuevo para el ACTON. ahora el arduino controla que pase q-switch
%   o no, así que lo voy a reprogramar '0'=no deja pasar, defecto;
%   '7'=activa relé y deja pasar q-switch, despues del 7º pulso
%       14dci2017 hay que subir NpulsesClean a 10 porque el ruido no se
%       estabiliza hasta el 6
%   v37 remove PIMAX/ACTON local and add code for the Avantes spectrometer.
%       As the PIMAX iCCD could be needed for some experiments, I'll set three global variables: USE_LOCAL_PIMAX, USE_AVANTES, USE_OCEANOPTICS
%   v38 shows the spectrum in real time to align the fiber
%   v39 remote capture of spectra with PIMAX3, needs the script
%   "remote_pimax.m" running in the host machine
%   v40 visualize remote pimax spectrum, try to reduce memory usage
%   v42 dual-delay: two delays for alternate capture for CFLibs NOT YET
%   IMPLEMENTED
%   v43 click&shot: option to choose the XY point of the next spatial point
%   v44 if USE_DUAL_PULSE_LASER==1 shutter and qswitching is controlled through serial
%   v45 add support for the STM32-based pulse generator
%   v46 the wavelength stitching of channels has changed
%   v47 the paramName value is parsed to perform real changes in every loop
%   v48 the upgrade to Windows10 changed the channel ordering, now master channel is #1
%   v49 trying to get "posX" loop names to go to absolute positions, besides changing deltaX
%   v50 fixed non-square pixels problem with manual_tracking at angles close to 45º.
%       We have found a variable qswitch delay approx. 400.800ns depending on the pump energy, I copy a polynomial fit  polyval([   -2.014032  107.14131 -1911.157624 11883.131128   ],pump)
%   v51 for new Amazon camera  (8mpx, ID=1)
%   v52 control of the absolute Z position
%   v53 extraCommands & new SCANLINE option=3
%   v54 add peak height & total intensity monitoring of spectra in preview window
%   v55 check brightness at the laser spot point to avoid measurement   outside the sample surface\
%   v57 recalibrated lambdas 7dic2021, the lambda vector is diferent
%   v58 store spectra8raw if needed, and the preview intensity bars are better suited to detect saturation
%       rules at 1mm & 0.1mm;  motor command "20 setjoybspeed" speeds up joystick when pressed
%   v59 11dic2022 remove calls to read position of motors just in case is not a memory problem
%   v60 17dic2022 handle sound with remote desktop, print X,Y coordinates, stitching, polygon
%   v61 22dic2022 major rework of movement: no initial reset to 0,0, motor position is read just once
%       added opengl('save','software') to get rid of the memory leak
%   v62 26dic2022 show sitched image in manual_tracking, LinePath implemented, new option NO_SHOW_IMAGE para evitar el memory leak
%       18ene2023 exit in case the Avantes spectrometer hangs up. fixed bug in starting_coords() when returning back from stitching
%   v63 new preview window with a new way of doing stitching 
%   v64 more friendly: some variables can be changed in the preview window, accel limit to 50mm/s2, new stitching with betterImage(TM), when loading stitched images, polygons and linepaths are effectively imported (shown in green)
%       25feb2024 version to process alumininum plates. extracommand is used to change parameters on-the-fly
%   v65 skipFunction = '' (do nothing), 'roi' (inROI as before) or 'anything' will call function @anything

opengl('save','software')
javaaddpath('C:\Program Files\Ocean Optics\OmniDriver\OOI_HOME\OmniDriver.jar'); % here because it CLEARS all the global vars
%global finish;
%global activar_qswitch;
%global pausing;
%global skipping;
%global lets_go;

global folderPath;
global experimentPath;
global datePath;
global filename;
global qswitch;

global vid; % camera objects are global and used by 'captureFrameGuppy.m'
global vid2; % webcam to look at pulse energy meter
global src;
global CAMERA; % camera number 1=guppy 2=amazon
global ID_CAMERA_SAMPLE; % if webcam plug in, sample camera ID is 2 % new v51 new camera from amazon 8mpx
global ID_CAMERA_PEM; % webcam if plugged

global pixel_shotX;  % needed by minimizer (image-based autofocus)
global pixel_shotY;
global focusWidth;
global focusHeight;
global focusLastZ;
global xyz;

dbstop if error; % this creates a breakpoint on errors and execution can be resumed with dbcont
disp('Starting...');
sound0=sin(0:0.1:500);sound1=sin(0:0.2:500);sound2=sin(0:0.4:1000);sound3=sin(0:1:2000);
pause('on');

%% CONFIGURATION THAT TIPICALLY IS NOT CHANGING IN EVERY EXPERIMENT
Npixeles_pimax = 1024;
Npixeles_hr2000 = 2048;
Npixeles_avantes = 15049; % 7dic2021
ValSaturation = 65000;   % expected maximum value without saturation
ValMinimun = 2000;
MAX_COLLECT_TRIES = 100;   % if no valid spectrum is collected, stop at this value
KEYPRESS_TO_CONTINUE = 0;  % if true, a key must be pressed to continue the next loop, it allows to prepare the laser again
SLEEP_BETWEEN_CAPTURES = 0.5; % pause to clean, think ... v55 We believe that some time is needed for the motors to stop, otherwise, positioning errors accumulate
MOTORS = 1;   % open or not COM: ports for motor control
BACKGROUND = 0; % capture or not the background signal, not very useful as a representative sample of the electronic background
RSD_PROCESSING = 1; %perform RSD online calculation
CAMERA = 2; % image capturing and processing. 1=guppy, 2=amazon %23dic2022 webcam removed to check usb hub power consumption
ID_CAMERA_SAMPLE = 2; % if webcam plug in, sample camera ID is 2
ID_CAMERA_PEM = 1; % webcam if plugged
ONLY_1PHOTO = 1; % only take a photo for param1Index=1 if Nparam2>1
TAKE_WEBCAM_PHOTO = 0; % 1= take webcam photo AFTER every experiment
SOUND = 0; %record audio
SOUND_STOP = 0;  % monitorizado del sonido y pausa cuando perfora
STD_SONIDO_STOP = 2.0; % valor umbral para parar, voya probar con el doble (x2), con porcus salía bien con 1.3 (30%)
PREVIEW = 1; % shows the video preview with the possibility of moving the target to a specific position
SPLIT_PARAM2 = 0; % if 1 , spectra() has only 1 param2 value, and after each param1 iteration the spectra is saved in an external file to save memory, v55 if 0, no split no matter how big are arrays
TWITTER = 0; % send tweets, 1=only when finishing, 2=when next direction for manual_tracking is needed, ...
ESCALA = 1; % draw 1mm ticks over the images
ORDENADOR_NUEVO=1;
USE_ARDUINO_QSWITCH = 0; % controls the qswitch signal using the arduino board, '0'=cut qsw signal, '1' let is pass, '7'=let it pass after 700ms
USE_SHUTTER_CTRL= 0; % if 1 the laser external shutter (thorlabs SC10 controller) is controlled by the program: close at startup, open and close around each experiment
USE_OCEANOPTICS=0; %=1 HR2000+ in hardware external trigger mode
USE_LOCAL_PIMAX=0; % it only has sense in the old computer where the PIMAX is installed
USE_AVANTES=1; % Avantes USB2 spectrometer (8 channels)
MULTI_WINDOW = 0; % realiza varias medidas a las lambdas indicadas en MW_centralwl. Guarda los espectros en MW_lambdas[:,num_ventanas] y MW_spectra(:,npulses,idx1,idx2,nw)
USE_REMOTE_PIMAX=0;   % simultaneous measurment with PIMAX3 in a remote computer
USE_DUAL_PULSE_LASER=1; % open COM12 to control laser lamp, shutter, qswitch and delays
USE_PULSE_GENERATOR=1; % uses the STM32-based pulse generator
USE_LEDLIGHT=1; % LED light controller with arduino in COM6, new v60 dic2022
CALIBRATED_PIXELS_MM=3;  % 1=lente 75mm, 2=objetivo x5, 3=objetivo x10
NO_SHOW_IMAGE=1; % no muestra la imagen de h_focus mientras mide, siempre quedan las fotos guardadas

%% CONFIGURATION THAT TIPICALLY CHANGES IN EVERY EXPERIMENT
FOLDER_NAME = 'Arqueo'; % folder: arqueo - aceros - carne - patatas - o whatever
EXPERIMENT_NAME =  'patellas_JAS';  % folder for this particular measurement 
SUBEXPERIMENT_NAME = '_nivel17-3';  % folder for this particular measurement 
EXPERIMENT_DESCRIPTION ='Medidas de patella vulgata 17.3  para el paper VS JAS';  % for the txt file
AUTOFOCUS_METHOD = 'image';  % methods:  'image' ; 'laserspot'
AUTOFOCUS = 0;  % perform autofocus between experiments. the initial position of the light spot is the reference, only de Z axis is changed.
% if (AUTOFOCUS)        h=warndlg('¡OJO! autofocus activo');     end;
DEFOCUS = 0;  %  if !=0 the sample is shifted in Z and back again after the measurement positive=up in mm
MANUAL_TRACKING = 0; % 1=after every experiment, show the video image and ask for a line to continue the scanning; 2=move just to that point
MANTRACK_RADIAL = 0; % if 1 the X movement is with angle motor, if 0, it is the X linear motor
DeltaT_ANGLE = 90; % degrees for the deltaT mode, 90 is transverse, for typical 2D mode
SCANLINE = 0; % =1:while shooting, the motor moves slowly along a line of length scanLineLength and angle scanLineAngle at the precise speed CENTERED AT THE INITIAL POINT; =2 the user draw the line and length or angle parameters are not valid. new parameter scanLineShrink change length in every iteration, should be 1 if not used; =3 same as 1 but movement start in the initial point, scanline is NOT centered
scanLineLength = 50; %% mm, in scanline=2 mode the length and angle is calculated from user rubberband input
scanLineAngle = 180; %% degrees 0=moves sample to the right, laser to the left; 180=move sample to the right
scanLineShrink = 1.0;  %% in each iteration, new length is calculated as previous length x shrink
AUTOCOLLECT = 0 ; % capture the numer of spectrum but one by one, ensuring all have peaks (this is useful for very low laser intensities close to the plasma formation threshold)
AUTOGAIN = 0; % if 1, if a spectrum is saturated, the actual gain value is multiplied by 0.9 and the new value recorded in experiment.txt
MEAN_INTENSITY = 8000; % mean intensity of the monitored lines
MAX_GAIN = 50;
USE_RECORDED_PATH = 0; % use a previously recorded path
ONLY_RECORDING_PATH = 0;  % if 1, laser is not shot and only the path is recorded and can be automated in manual_tracking mode
recorded_path_filename = 'D:\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\21122022-114045_LAN176_seccion-izquierda_lado-anterior_PATH'; % <<<<<<<<<<<<<<<<<<<<<   CHANGE PATH HERE WITHOUT THE FILENAME, JUST COPY&PASTE
CHECK_SURFACE = 0; % if 1, laser is not shot if brightness at the focal point is bellow a threshold; if =2, DOES NOT END PARAM1 IF A DARK POINT IS FOUND (that is, if a dark zone is found after a bright path, it goes on until the end
checkSurface_black = 60; % this is the TRHESHOLD  brightness , below that value, it is assumed there is no surface at the focal point
STORE_SPECTRA8RAW=0; % if true, the complete CCD spectra are stored in "spectra8raw" variable, a warning is issued if RAM usage is too high
MinXmmStitch=1; % area to be stitched
MaxXmmStitch=2;
MinYmmStitch=1;
MaxYmmStitch=2;
dxStitch=0.6; % the x,y witdh,height of each tile, better no to change, as it is limited by distortion of camera image
dyStitch=0.4;
havePolygon = 0; % ROI captured. if 1, then the file below is loaded before measuring
filenamePolygon = 'D:\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\';
haveLinePath = 0; % equivalent to measuringPath
filenameLinePath = 'D:\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\';
haveImageStitched = 0; % load a previously saved figure with the stitching
filenameImageStitched = 'D:\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\Lapas_Langre2012_Repulidas\';

%parámetros de medida
LENTEFIBRA = 'Fibra desnuda 1mm';%2 lentes nuevas + 1mm + roja Avantes'; % 'PIMAX3: azul600um sin lente, HR2000: con lente'
BOMBEO= 'variable pumpEnergy 16jul2024: sigue igual, 18J = 21mJ antes de espejo ';
FOCO = 'en el foco (enfoque de la imagen) con objetivo x10';
BICHO= '';
MUESTRA='';
RECOGIDA = '';
TIPOMUESTRA = '';
NUMPASADA = '';

%capture
Npulses = 15; %  laser pulses per point (per experiment)
NpulsesClean = 4; % is in fact for noise capture, laser is not activated!! if 3 or more, these N first pulses are discarded (laser qswitch not triggered). 7= 5 pulses out guaranteed
frequencyPulses = 10;  % Hertz
gainIntensifier = 10; % gain only affects PIMAX3
UseAlternateDelay= 0; % 1=alternate delays in consecutive shots, only with STM32 custom pulse generator; 2=alternate delays in consecutive spatial ponts, done by hand   NO YET IMPLEMENTED
AlternateDelay = 0; % [ns] of alternate delay (only if USE_PULSE_GEN==1)
delayCaptureForReal = 1000;  % [ns]  REAL delay between laser pulse and spectrometer capture window, calculated for every setup
pumpEnergy = 18; % 16jul2024: sigue igual, 18J = 21mJ antes de espejo 26feb2024:   dualpulse 15.8J = 10mJ antes de espejo;  17.8J=20mJ ; 23.5J=50mJ (demasiado) era 14.9 con la lente habitual; 9feb2023: 17J=28mJ antes del pinhole
DualPulse = 1; % activates second channel on DP lotis laser
DualPulseDelay = 1000; % [ns]  channel 2 shots this ns BEFORE the channel 1 OJO!! reducido porque salía más intensidad para retrasos cortos
gateWidth = 20; % en nanosegundos/microsegundos OJO!! he comprobado que poniendo "20" lo toma como 20us, pero poniendo "500" lo toma como 500ns, poniendo "0.5" lo toma como 0.5 nanosegundos , "1000" son 1000 nano(!)
gains_set = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]; % gains for every pellet
avantesIntegrationTimems = 1.05; % integration time for avantes

%lineMonitoring_name = {'CaI@227'; 'CaI@239';'CaI@558';'CaII@315';'MgI@517';'MgII@280'}; % new v54 those lines heights are monitored in the "disparar" preview window.
%lineMonitoring_lambda = [ 227.501 , 239.78, 558.162, 315.85, 517.234, 280.228]; % OJO!! lambdas descalibradas al 4nov2021
%lineMonitoring_min = [ 1150, 1400, 20000, 9000, 6000, 2900];
%lineMonitoring_max = [ 1250, 1600, 33000, 15000, 11000, 4800];
%3dic2021 ALL lines
lineMonitoring_name = {'CaI@227'; 'CaI@239';'CaI@558';'CaII@315';'CaII@370';'CaII@373';'MgI@382';'MgI@516';'MgI@517';'MgI@518';'MgII@279';'MgII@280'}; % new v54 those lines heights are monitored in the "disparar" preview window.
%lineMonitoring_lambda = [ 656.279, 227.501 , 239.78, 558.162, 315.85, 370.603, 373.69, 382.936, 516.732, 517.2684, 518.36, 279.55, 280.27]; % originales, usar estas cuando esté recalibrado
lineMonitoring_lambda = [ 227.546;  239.856; 558.197; 315.89; 370.603; 373.690; 382.9355; 516.7322; 517.2684; 518.360; 279.55; 280.2704 ]; % lambdas LA-ICP with CALIBRATED lambdas
lineMonitoring_min = [   900,    850,      20000,   4000,   10000,   12000,  1600,    3000,    7000,     12000,  4000,   2000];
lineMonitoring_mean = [  972,    950,     29312,    10705,  26625,   35063,  2034,    4451,    11885,    19093,  7023,   3595];
lineMonitoring_max = [  1100,    1050,    45000,    17000,  42000,   52000,  2500,    7000,    18000,    27000,  10000,  5000];



%%%%%%%%%%%%%%%%%%%%%  names with specific meaning     new! v47
%   'deltaT' transverse movement (2D scanning) for eachpoint, 'DeltaT_ANGLE' is the angle
%   'deltaR' moves along a path, it autom. changes variable deltaR=param1Step and set MANUAL_TRACKING=1
%   'click&shot' set MANUAL_TRACKING=2
%   'DualPulseDelay_ns'  the delay between laser pulses (ns) is changed, set dualpulse=1
%   'delayCaptureForReal_ns' the espectrometer capture delay (ns) is changed within the loop. changes are written to avantes and pimax spectrometer
%   'posX' set the variable deltaX=param1step AND go to absolute position X,Y = posX,posY based on idx values, to perform a spatial scanning
%   'posY' set the variable deltaY=param1step AND go to absolute position based on idx values, to perform a spatial scanning
%   'posZ' set the variable deltaZ=param1step AND go to ABSOLUTE position, to perform a depth profiling
%   'gain' set the gain for PIMAX camera
%    new v53  extracommands perform additional actions right before starting the laser shot burst in each spatial point
%       extraCommands = 'fprintf(xyz,'' 1 0 r\n'')'; % move sample 1mm to the left, laser to the right
%    medida 2D ¡OJO! si avanza en negativo, poner step negativo
extraCommands = '';
%extraCommands = 'pause;'; % do nothing
%extraCommands = 'fprintf(xyz,''0 0.4 0 r\n'');' % move +1mm x axis relative
skipFunction = 'roi'; % skipFunction = '' (do nothing), 'roi' (inROI as before) or 'anything' will call function @anything

param1Name = 'posX'; param1From = 0; param1To = 2; param1Step = 0.05;  % grid scanning
param2Name = 'none'; param2From = 0; param2To = 0; param2Step = 0.05;

%barrido en Z para ver el foco y diametro de spot
%param1Name = 'posZ'; param1From = -0.5; param1To = 0.5; param1Step = 0.1;  % grid scanning
%param2Name = 'posY'; param2From = 0; param2To = 0; param2Step = 0.3;

%param1Name = 'posX'; param1From = 0; param1To = 0.9; param1Step = 0.05;  % grid scanning
%param2Name = 'posY'; param2From = 0; param2To = 0.4; param2Step = 0.05;

%param1Name = 'DualPulseDelay_ns'; param1From = 500; param1To = 5000; param1Step = 500;  % linear path
%param2Name = 'none'; param2From = 0; param2To = 0; param2Step = 1;

%param1Name = 'cordon'; param1From = 1; param1To = 1; param1Step = 1;  %
%param2Name = 'none'; param2From = 1; param2To = 1; param2Step = 1;  % none

%param1Name = 'punto'; param1From = 1; param1To = 5; param1Step = 1;  %
%param2Name = 'none'; param2From = 1; param2To = 1; param2Step = 1;  % none

%param1Name = 'posZ'; param1From = 0; param1To = 5; param1Step = 0.2; %recorrido en Z, buscar punto focal
%param2Name = 'none'; param2From = 1; param2To = 1; param2Step = 1;  % none

%param1Name = 'posX'; param1From = 0.3; param1To = 0.5; param1Step = 0.2;  % grid scanning
%param2Name = 'posY'; param2From = 0; param2To = -3.8; param2Step = 0.1;

%param1Name = 'deltaT'; param1From = -0.2 ; param1To = 0.2; param1Step = 0.05;  %  9pts Npoints2D = (param1To-param1From)/param1Step + 1;
%param2Name = 'deltaR'; param2From = 0; param2To = 20; param2Step = 0.05;  %

%param1Name = 'pulsos_conZ'; param1From = 0 ; param1To = 2000; param1Step = 50;  %  poniendo 'pulsos' no se mueve (a propósito)
%param2Name = 'none'; param2From = 1; param2To = 1; param2Step = 1;  % none

%param1Name = 'delay'; param1From = 800 ; param1To = 1000; param1Step = 50;  %  poniendo 'pulsos' no se mueve (a propósito)
%param2Name = 'none'; param2From = 1; param2To = 1; param2Step = 1;  % none

%param1Name = 'deltaX'; param1From = 0 ; param1To = 2; param1Step = 0.2;  %  poniendo 'pulsos' no se mueve (a propósito)
%param2Name = 'none'; param2From = 1; param2To = 1; param2Step = 1;  % none


%movement  EY!! OJO!! coordinates changed from physical to image
deltaX = 0.0;   % relative displacement in X axis after each point IMAGE COORDINATES!!! OJO!! 4setp2014: he cambiado los signos en el sprintf al motor porque estaban al revés cuando no había manual_tracking
deltaY = 0.0;  %% ojo eje Y tambnién invertido en no manual_tracking -> cambiado de signo a ver si se arregla
deltaZ = 0.0;   % positive is upwards
deltaP = 0.0; % mm of displacement using the Angle motor; POSITIVE->upper surface to the real right
deltaR = 0.05; % used by MANUAL_TRACKING feature, for BOTH linear and radial movement
DIAMETER = 16.2;   % diameter in mm of the rotating sample, it is the diameter at the START of an experiment
DIAMETER_END= 14.8;  % the diameter of the sample at the END of the path, the angle changes linearly down to this diameter
everyNpoint = 9999; % AFTER nth point the distance deltaR and deltaP is change by everyIncDec, if not used, set to 9999
everyIncDec = 0.25; % the distance every nth point is incremented by this value, in mm (only affect deltaR and deltaP)


%lines selected for saturation detection AND ratio. If satlines[] is not empty, their
%pixels are compared against ValSaturation to discard spectra
%satlines = [ 277.9831 279.5528 280.2704 285.2127 300.686 299.496 299.731 299.964];  % líneas interesantes espectro 1200lpp @ 290nm original
%satlines = [ 279.49 280.27 285.21  315.89 370.603 373.69 383.83 487.813 504.16 518.88  ];  % líneas interesantes espectro 150lpp @ 420nm
%satlines = [ 383.83 487.813 504.16 518.88  ];  % líneas interesantes espectro 150lpp @ 420nm
%satlines = [ 570.35 585.75  ];  % DESCALIBRADO líneas interesantes en la ventana 1800lpp 578nm
%satlines = []; %is empty, the global maximum is compared with ValSaturatio
%satlines = [ 535.2 570.8  ];  % una de calcio y otra de magnesio
%satlines = [ 383.83 395.705 ];  % Mg y Sr
%satlines = [ 571.1 585.75  ];  % líneas interesantes en la ventana 1800lpp 578nm
%satlines = [ 383.83 370.603 ];  % buscando líneas buenas para los porcus
%satlines = [ 285.15 300.64 ];  % líneas interesantes en la ventana 150lpp 416nm
%satlines = [ 407.7 428.3 ];  % líneas interesantes en la ventana 1200pp 418nm
satlines = [ 460.7 670.7 ];  % líneas interesantes en la ventana 150lpp 416nm  [ 285.2 300.7 ];
%satlines = [ 407.7 428.3 ];  % estroncio calcio en la ventana 1200pp 397nm

LAMBDA_MAGNESIO = [285.2127]; %MgI
LAMBDA_CALCIO = [300.686];% CaI
background4lambdas_mg = [ 282 284.6 285.7 289.3 ];  %lambdas_MgI en el pico de 285.2
background4lambdas_ca = [ 295.4 299 301.3 305.3  ]; %ñambdas para lambda 300.686
LAMBDA_RESINA = 288.1; % silicio
background4lambdas_resina = [ 287 287.8  288.3 289.3  ]; %lambdas para resina en 288.1

lineMonHeights = zeros(numel(lineMonitoring_lambda),1); % array of values to show
lineMonLabels = cell(numel(lineMonitoring_lambda),'');
lineMonString = strings;

for i=1:numel(lineMonitoring_lambda)
    lineMonLabels(i) = lineMonitoring_name(i)  ;
end

%MW_centralwl = [ 290, 322 , 354, 386, 418, 449, 480, 510, 540, 569];  % 1200lpp, a partir de 480 la eficiencia es muy baja
%MW_centralwl = [ 280, 312 , 344, 376, 408, 439 ];  % 1200lpp, a partir de 480 la eficiencia es muy baja
MW_centralwl = [ 290, 322 , 354, 386, 418, 449 ];  % 1200lpp, a partir de 480 la eficiencia es muy baja



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NO HAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QUE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMBIAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NADA MÁS


blankImage = zeros(960,1280,3,'uint8'); % v05 no camera = blank image
if isequal(skipFunction,'')
    skipFunctionHandler = 0;
else    
    skipFunctionHandler = str2func(skipFunction);
end
StitchedIsLoadedNotCreated=0;
if haveImageStitched
    try
        h_stitch=openfig(filenameImageStitched);
        disp(['Loaded ' filenameImageStitched]);
        StitchedIsLoadedNotCreated=1;
        p=findobj(h_stitch,'Type','images.roi.polygon'); % check if polygon is already drawn
        if ~isempty(p)
            roiPol = drawpolygon('Position',p(1).Position);
            delete(p); % new v17
            set(roiPol,'Color','green');
        end
        p=findobj(h_stitch,'Type','images.roi.polyline'); % check if linepath is already drawn
        if ~isempty(p)
            roiLine = drawpolyline('Position',p(1).Position);
            delete(p); % new v17
            set(roiLine,'Color','green');
        end
    catch
        error('ERROR!! the filename with the stitchedImage cannot be read, check name');
    end
    h_stitch.Name = 'Loaded stitched image';
    %buttons are not used but needed if disabled by manual_tracking
    esc_button_stitch = uicontrol('Style', 'PushButton', 'String', 'Finish', 'Callback', 'finishStitch=1;', 'Position',[10 10 50 20], 'ForegroundColor','red','Enable','off');
end

if havePolygon & haveImageStitched % if was set to 1 before, then we should load the roi file now
    try
       load(filenamePolygon);
       disp(['Loaded ' filenamePolygon]);
       StitchedIsLoadedNotCreated=1;
       roiPol = drawpolygon('Position',roiPol.Position);
    catch
        error('ERROR!! the filename with the Polygon cannot be read, check name');
    end
end
if haveLinePath % if was set to 1 before, then we should load the roi file now
    try
       load(filenameLinePath);
       disp(['Loaded ' filenameLinePath]);
       StitchedIsLoadedNotCreated=1;
       roiLine = drawpolyline('Position',roiLine.Position);
    catch
        error('ERROR!! the filename with the LinePath cannot be read, check name');
    end
end


%capture delays calculated for the different pulses
DualPulse_Channel1_delay = 160000; % fixed for optimum lamp charge
DualPulse_Channel2_delay = 160000 - DualPulseDelay; % fixed for optimum lamp charge
% 20ago2019 with external lamp (pulsegen) and internal qswitches (160us) the laser shot is at 160900 ns (21ago2019)
% 2sept2019 with external lamp and external qswitches 160us -> pulso qswitch sale en 161700ns, the laser shot is
% exactly at 162000ns
%OJO!!! todos los tiempos están mal , hay que recalibrar el stm32, lo dejo
%apañado el 3sept2019 con lso pulsos verificados, con 159000 genera el
%pulso del espectrómetro a 160700, +1300 = 162000 que es cuando dispara el
%láser
% el 6setp2019 verificado laser shot at 160400ns with external qswitch,
% el 21ene2021 verificado laser shot at 160800ns with external qswitch, its variable
% avantes capture delay = 1300ns
delayCaptureAvantes = 160000 - 1300 + delayCaptureForReal;  % [ns]
delayCaptureAvantes = delayCaptureAvantes + round(polyval([   -2.014032  107.14131 -1911.157624 11883.131128 ],pumpEnergy)); % polynomial fit of qswitch delay

%capture for REMOTE PIMAX only
remote_pimax_npulses = Npulses; % +6; % 2setp2019 there are six pulses with no avantes sync, but number of pulses should be Npulse
remote_pimax_gainIntensifier = gainIntensifier;
remote_pimax_delayCapture = delayCaptureAvantes + 1300; % el avantes tarda realmente 1300ns en capturar desde el trigger
remote_pimax_gateWidth = gateWidth; %
absoluteMovement = 0; % if 1 (using posX, posY or posZ) the movements are ABSOLUTE according to loop parameters, and performed BEFORE measurement

% first-time changes of variables or whatever for specific names in the loops
if ( strcmp(param1Name,'deltaR') == 1 || strcmp(param2Name,'deltaR') == 1)
    deltaR = param1Step; % set the step, in mm
    MANUAL_TRACKING = 1; % activate the manual tracking mode
end
if ( strcmp(param1Name,'click&shot') == 1 || strcmp(param2Name,'click&shot') == 1)
    MANUAL_TRACKING = 2; % activate the manual tracking 'click & shot'
end
if ( strcmp(param1Name,'DualPulseDelay_ns') == 1 || strcmp(param2Name,'DualPulseDelay_ns') == 1)
    DualPulse = 1; % set dualpulse
end
if ( strcmp(param1Name,'posX') == 1)
    deltaX = param1Step; % set de X movement to the step of the loop
    absoluteMovement = 1;
end
if ( strcmp(param2Name,'posX') == 1)
    deltaX = param2Step; % set de X movement to the step of the loop
    absoluteMovement = 1;
end
if ( strcmp(param1Name,'posY') == 1)
    deltaY = param1Step; % set de Y movement to the step of the loop
    absoluteMovement = 1;
end
if ( strcmp(param2Name,'posY') == 1)
    deltaY = param2Step; % set de Y movement to the step of the loop
    absoluteMovement = 1;
end
if ( strcmp(param1Name,'posZ') == 1)
    deltaZ = param1Step; % set de Z movement to the step of the loop
    absoluteMovement = 1;
end
if ( strcmp(param2Name,'posZ') == 1)
    deltaZ = param2Step; % set de Z movement to the step of the loop
    absoluteMovement = 1;
end
if absoluteMovement==1
     MANUAL_TRACKING = 0; % just in case, I tend to forget
end
%to activate 2D scanning, param1 (or param2?) should be deltaT (T for transverse)
if ( strcmp(param1Name,'deltaT') == 1)
    Npoints2D = round((param1To-param1From)/param1Step + 1);  % number of total points of the transverse line
    L2D = (param1To-param1From);  % milimeters of the scanning line
elseif ( strcmp(param2Name,'deltaT') == 2)
    Npoints2D = round((param2To-param2From)/param2Step + 1);
    L2D = (param2To-param2From);
else
    Npoints2D=1 ; % not 2D scanning
end


if ORDENADOR_NUEVO  % mapeado de puertos en el nuevo
    COM_XYZ='COM4'; % es el puerto COM4 en la placa nueva (oct2017), era COM2 en giflibs, ha cambiado a COM10, y en algún momento de 2020 cambió a com4
    COM_THETA='COM1'; % puerto interno de la placa conectado al pic de luis (theta)
    COM_SHUTTER='COM8'; % shutter from thorlabs, prolific USB-serial adaptor (Adolfo)
    COM_QSWITCH='COM7';
    COM_DUAL_PULSE='COM12';
    COM_PULSE_GENERATOR='COM13';
    COM_LEDLIGHT='COM6'; % new v60 LED controller
else % mapeado de puertos en el viejo
    COM_XYZ='COM2'; % COM2 en giflibs antiguo
    COM_THETA='COM1'; % puerto interno de la placa conectado al pic de luis (theta)='COM16'; % shutter from thorlabs, prolific USB-serial adaptor (Adolfo)
    COM_QSWITCH='COM5';  % "arduino uno"
end

remote_pimax_ip = '193.144.200.68'; % giflibs.teisa.unican.es

%twitter
if (TWITTER)
    credentials = struct('ConsumerKey' , 'xxx', 'ConsumerSecret' ,  'xxx', 'AccessToken' , 'xxx' , 'AccessTokenSecret' , 'xxx');
    try
        tw = twitty(credentials);
    catch
        TWITTER=0;
    end
end

%computer
folderPath =  ['D:\UNICAN\Proyecto deepRAMP - LIBS\' FOLDER_NAME '\'];
experimentPath = [EXPERIMENT_NAME,'\'];
backgroundFile = 'background.spe';
dataFile = 'matlabData.mat';
if (USE_RECORDED_PATH ==1 && strcmp(recorded_path_filename,'')) % a specific path is not specified
    recorded_path_filename = strcat(folderPath,experimentPath,datePath,'recorded_path.mat');
else
    recorded_path_filename = strcat(recorded_path_filename ,'\recorded_path.mat');
end
if (ONLY_RECORDING_PATH)
    datePath = strcat(datestr(now,'ddmmyyyy-HHMMSS'),SUBEXPERIMENT_NAME,'_PATH\');
else
    datePath = strcat(datestr(now,'ddmmyyyy-HHMMSS'),SUBEXPERIMENT_NAME,'\');
end
setupFile = 'setup.exs';

%autofocus
ref_pixel = [0 0];
act_pixel = [0 0];
%calibrated on 11jun2014 with blue laser over p. lineatus surface
calibrated_pixels_mm_Z = 376; % calibrated pixels per Z mm displacement; positive displacement = positive pixel increase (camera china-luis 26feb2015)
%calibrated_pixels_mm_Z = 182; % 12jul2018 como el láser está en ángulo,
%depende del ángulo en sí de la muestra
focusWidth = 100;
focusHeight = 100;
focusMethod = 'GRAS'; % gradient

%pixel_shotX = 659;  % coordinates of the laser shot over the image, now v33 read from .mat file, created by the "calibrate laser-spot" code
%pixel_shotY = 477;
load('last_laserspot_coordinates.mat'); %new v33

%calibrated_pixels_mm_XY = 82; % 82 pixels = 1 mm in XY plane. GUPPY

if CALIBRATED_PIXELS_MM==1     %calibración para la lente de 75mm y cámara Amazon (a 2nov2022)
    calibrated_pixels_mm_XY_X = 1259/7; %  11ene2021: parece que la resolución 1280x1024 que tenia no tiene los píxeles cuadrados, lo desdoblo en XY, 75mm lens
    calibrated_pixels_mm_XY_Y = 898/5; %   9feb2021: new amazon camera, square pixels, thanks god, 675mm lens
elseif CALIBRATED_PIXELS_MM==2
    calibrated_pixels_mm_XY_X = 325.7; %  15jun2021: new amazon camera, x5 NIR objective
    calibrated_pixels_mm_XY_Y = 325.7; %
elseif CALIBRATED_PIXELS_MM==3      %calibración para objetivo x10 NIR
    calibrated_pixels_mm_XY_X = 658; %  2nov2022: new amazon camera, x10 NIR objective
    calibrated_pixels_mm_XY_Y = 658; %  trying 668 instead of 658 17jul2024, it does not change anything in the stitching maps
else
    error('ERROR: opción de calibración pixels/mm que no existe');
end

mantrack_next_points = 0; % number of points that can shot automatically
mantrack_deltaX = 0; % increment for the next point in a manual tracking sequence
mantrack_deltaY = 0;
recorded_path = zeros(1,5); % record the trayectory for manual or automatic tracking: X, Y, Z, autofocus, deltaP; it size increments with each point because the number of points is not known
finish = 0;
activar_qswitch = 0;
pausing = 0;
skipping = 0;
energia_limit = 100; % for the spectral energy bar
total_ws_steps = 0; % keep a count of forward (W) and backwards (S) stepper motor steps to get back to the origin
total_distance = 0.0; % keep track of total deltaR/deltaP distance
%recorded path
if (USE_RECORDED_PATH && ~haveLinePath) % new v62
    load (recorded_path_filename);
    old_recorded_path = recorded_path; % new recorded path is updated anyway, it should be the same
end
%audio
if (SOUND)
    try
        recorder = audiorecorder(44100, 16, 1);
    catch
        disp('WARNING!! sound recording deactivated - wont work using remote desktop');
        SOUND=0;
    end
end

%comprobar los puertos serie
%% v31 comprueba los puertos serie disponibles
serialInfo = instrhwinfo('serial'); % requires instrument control toolbox
if (MOTORS)
    puerto=strfind(serialInfo.SerialPorts,COM_THETA);
    if (isempty(cell2mat(puerto))) % no existe el puerto COM_THETA
        error('¡ATENCIÓN! el puerto COM_THETA del controlador del motor angular no está disponible.')
    end
    puerto=strfind(serialInfo.SerialPorts,COM_XYZ);
    if (isempty(cell2mat(puerto))) % no existe el puerto COM2
        error('¡ATENCIÓN! el puerto COM2 del posicionador XYZ no está disponible.')
    end
    
end
if (USE_SHUTTER_CTRL)
    puerto=strfind(serialInfo.SerialPorts,COM_SHUTTER);
    if (isempty(cell2mat(puerto))) % no existe el puerto COM6
        error('¡ATENCIÓN! el puerto COM6 del shutter no está disponible.')
    end
end
if (USE_ARDUINO_QSWITCH)
    puerto=strfind(serialInfo.SerialPorts,COM_QSWITCH);
    if (isempty(cell2mat(puerto))) % no existe el puerto COM5
        error('¡ATENCIÓN! el puerto COM5 del sincronizador del Q-switch (o sea, el arduino) no está disponible. Apretar botón 1 del hub USB.')
    end
end

if (USE_DUAL_PULSE_LASER)
    puerto=strfind(serialInfo.SerialPorts,COM_DUAL_PULSE);
    if (isempty(cell2mat(puerto))) % no existe el puerto COM5
        error('¡ATENCIÓN! el puerto COM12 del control del láser LOTIS dual-pulse no está disponible. ¿esta el láser apagado?')
    end
end

if (USE_PULSE_GENERATOR)
    puerto=strfind(serialInfo.SerialPorts,COM_PULSE_GENERATOR);
    if (isempty(cell2mat(puerto))) % no existe el puerto
        error('¡ATENCIÓN! el puerto COM13 del generador de pulsos (STM32) no está disponible');
    end
end

if (USE_LEDLIGHT)
    puerto=strfind(serialInfo.SerialPorts,COM_LEDLIGHT);
    if (isempty(cell2mat(puerto))) % no existe el puerto COM6
        error('¡ATENCIÓN! el puerto COM6 del control de iluminación no está disponible.')
    end
end


%% WixX32 objects
if (USE_LOCAL_PIMAX)
    objCal = actxserver('WinX32.CalibObj');
    objWin = actxserver('WinX32.DocWindows');
    objDocs = actxserver('WinX32.DocFiles');    % It seems needed to access the calibration of the first DOC
    objDoc = actxserver('WinX32.DocFile');    %
    objExp = actxserver('WinX32.ExpSetup');
    objSpe = actxserver('WinX32.SpectroObj');
    objPul = actxserver('WinX32.PITG');
    objSpes = actxserver('WinX32.SpectroObjMgr'); % se necesita para mover el grating
end

%%HR2000
if (USE_OCEANOPTICS==1)
    %% probar high-speed-acq
    %% Connect the spectrometer
    %   try
    wrapper = com.oceanoptics.omnidriver.api.wrapper.Wrapper();
    NoOfDevices = wrapper.openAllSpectromete
    if (NoOfDevices >1)
        h=warndlg('¡ATENCIÓN! hay más de un espectrómetro HR2000+ conectado y se utilizará solo el número 0');
    end
    %   catch
    %      error 'Se ha activado HR2000=1 y no se ha podido inicializar el wraper ¿algo de java?¿no se ha instalado?'
    %   end
    
    %% configuración por defecto
    wrapper.setTimeout(0,5000); % timeout for trigger mode 4
    wrapper.setIntegrationTime(0,1000); %1000 = 1ms
    wrapper.setExternalTriggerMode(0,4); % hardware external trigger, mode #4 after firmware revision
    wrapper.setScansToAverage(0,1);
    hr2000_lambdas = wrapper.getWavelengths(0); % el HR2000+ tiene 2048 lambdas
end

% here there was the Avantes initialization code, moved down so
% spectrometers are free to use by its software
% v38 moved up again so the spectrum can be seen in the preview window, the
% avantes software is no longer needed
% v48 windows upgrade changed the channel numbering, now master is #1
if (USE_AVANTES)
    global h1 h2 h3 h4 h5  h6 h7 h8 S1
    
    nport_avantes=spectrometer('init',0); % number of spectrometers
    if (nport_avantes < 8)
        error ('El espectrómetro avantes no es accesible ¿cerrar el software de Avantes?');
    end
    ID_avantes=spectrometer('getlist'); % data of each one
    
    h5=spectrometer('activate',ID_avantes(5)); % 400-490nm
    h6=spectrometer('activate',ID_avantes(6)); % 490-560nm
    h7=spectrometer('activate',ID_avantes(7)); % 560-680nm
    h8=spectrometer('activate',ID_avantes(8)); % 680-880nm
    h1=spectrometer('activate',ID_avantes(1)); % 180-260nm MASTER
    h2=spectrometer('activate',ID_avantes(2)); % 251- 316
    h3=spectrometer('activate',ID_avantes(3)); % 312 - 365
    h4=spectrometer('activate',ID_avantes(4)); % 360 - 403
    spectrometer('usehighres',h1,true);
    spectrometer('usehighres',h2,true);
    spectrometer('usehighres',h3,true);
    spectrometer('usehighres',h4,true);
    spectrometer('usehighres',h5,true);
    spectrometer('usehighres',h6,true);
    spectrometer('usehighres',h7,true);
    spectrometer('usehighres',h8,true);
    % measurement configuration
    S1.TriggerMode=1; % 0=software, 1=hardware
    S1.TriggerSource=0; % 0=synchronized 1=external ! it should be external with the cable connecting pin6 of DB26, but it does not work
    S1.TriggerSourceType=0; %0=edge 1=level
    S1.IntegrationTime=avantesIntegrationTimems; % hang-up with values below 2.0 ¿?  28mayo2018 no, 1.05 is the minimum
    S1.StartPixel=0;
    S1.StopPixel=2047;
    S1.IntegrationDelay=0;  % check if this is "laser delay"
    S1.CorDynDark=0;
    S1.Smoothing=0;
    S1.SaturationDetection=1;
    %myLambda1=spectrometer('getlambda',h1);    spectrometer('measconfig',h1,S1); % it is the master, last to measure
    myLambda2=spectrometer('getlambda',h2);    spectrometer('measconfig',h2,S1);
    myLambda3=spectrometer('getlambda',h3);    spectrometer('measconfig',h3,S1);
    myLambda4=spectrometer('getlambda',h4);    spectrometer('measconfig',h4,S1);
    myLambda5=spectrometer('getlambda',h5);    spectrometer('measconfig',h5,S1);
    myLambda6=spectrometer('getlambda',h6);    spectrometer('measconfig',h6,S1);
    myLambda7=spectrometer('getlambda',h7);    spectrometer('measconfig',h7,S1);
    myLambda8=spectrometer('getlambda',h8);    spectrometer('measconfig',h8,S1);
    % master channel should be configured as external?
    %S1.TriggerSource=1;% 0=synchronized 1=external it does not work
    myLambda1=spectrometer('getlambda',h1);    spectrometer('measconfig',h1,S1); % same configuration for all
    %new 7dic2021: dll calls do not get the new wavelengths coefficients,I'll do it
    px=0:2047;px=px';
    l1 = -3.146179E-6.*px.*px +4.401482E-2.*px + 1.785329E2;
    l2 = -1.643487E-10.*px.*px.*px -3.171126E-6.*px.*px + 3.916221E-2.*px + 2.513996E+2;
    l3= 6.156498E-10.*px.*px.*px -4.122094E-6.*px.*px + 3.400006E-2.*px + 3.123916E2;
    l4= -4.901761E-10.*px.*px.*px -3.092496E-6.*px.*px + 2.898256E-2.*px + 3.608379E2;
    l5= -2.034117E-10.*px.*px.*px -4.94692E-6.*px.*px + 5.679243E-2.*px + 3.986681E2;
    l6= -3.812122E-10.*px.*px.*px -5.009936E-6.*px.*px + 4.896741E-2.*px + 4.889888E2;
    l7= -2.451990E-10.*px.*px.*px -7.002403E-6.*px.*px + 7.403560E-2.*px + 5.602282E2;
    l8= -5.067776E-10.*px.*px.*px -8.812136E-6.*px.*px + 1.216518E-1.*px + 6.816069E2;
    avantes_lambdas = cat(1,l1(1:1921),l2(1:1867),l3(1:1682),l4(1:1694),l5(1:1930),l6(1:1859),l7(1:2048),l8(1:2048));
    %version < 7dic2021
    %avantes_lambdas = cat(1,myLambda1(1:1921),myLambda2(1:1867),myLambda3(1:1786),myLambda4(1:1755),myLambda5(17:1943),myLambda6(1:1859),myLambda7(1:2048),myLambda8(1:2048));
    %version 7dic2021
    %avantes_lambdas = cat(1,myLambda1(1:1921),myLambda2(1:1867),myLambda3(1:1682),myLambda4(1:1831),myLambda5(17:1930),myLambda6(1:1859),myLambda7(1:2048),myLambda8(1:2048));
end

%% REMoTE_PIMAX
if (USE_REMOTE_PIMAX==1)
    try
        % parameters of capture
        ahora1=tic;
        sock = msconnect(remote_pimax_ip,3000);
        comando = [ 'npulses=' num2str(remote_pimax_npulses)];
        ret=mssend(sock, comando, 2);ack = msrecv(sock,3); msclose(sock);
        ahora2=toc(ahora1);
        if (ahora2 > 3)
            disp('ERROR! Se ha especificado REMOTE_PIMAX=1 pero el ordenador remoto no es accesible');
            error('We better stop now.');
        end
        sock = msconnect(remote_pimax_ip,3000);
        comando =  ['gain=' num2str(remote_pimax_gainIntensifier)];
        ret=mssend(sock, comando, 2);ack = msrecv(sock,3); msclose(sock)
        sock = msconnect(remote_pimax_ip,3000);
        comando = [ 'delay=' num2str(remote_pimax_delayCapture)];
        ret=mssend(sock, comando, 2);ack = msrecv(sock,3); msclose(sock);
        sock = msconnect(remote_pimax_ip,3000);
        comando = [ 'width=' num2str(remote_pimax_gateWidth)];
        ret=mssend(sock, comando, 2);ack = msrecv(sock,3);msclose(sock);
        if (strcmp(ack,'ok')==0)
            warndlg('¡ATENCIÓN! no hay respuesta del servidor REMOTE_PIMAX');
        end
    catch
        disp('Error TCP/IP accediendo a la cámara PIMAX3 remota.');
    end
end
%% Start

if (MOTORS)
    %open serial port for motor Theta at COM_THETA: and XYZ (Corvus XYZ, positioner) at COM_XYZ:
    try
        %theta=serial(COM_THETA); theta.BaudRate=38400;theta.Terminator='CR';fopen(theta);if(theta.BytesAvailable>0);fread(theta,theta.BytesAvailable);end
        xyz=serial(COM_XYZ); xyz.BaudRate=57600;fopen(xyz);xyz.FlowControl='software';if(xyz.BytesAvailable>0);fread(xyz,xyz.BytesAvailable);end
        % example:  fprintf(xyz,' -1 0 0 r'); move -1mm x axis relative
        % if not properly closed, exec:  delete(instrfindall)
        % serialInfo = instrhwinfo('serial')
        % IMPORTANT! theta serial port only works with fwrite(theta,'w'); fprintf
        % is not working!
    catch
        delete(instrfindall);
        disp('Warning! serial ports not properly closed');
        %theta=serial(COM_THETA); theta.BaudRate=38400;theta.Terminator='CR';fopen(theta);if(theta.BytesAvailable>0);fread(theta,theta.BytesAvailable);end
        xyz=serial(COM_XYZ); xyz.BaudRate=57600;fopen(xyz);xyz.FlowControl='software';if(xyz.BytesAvailable>0);fread(xyz,xyz.BytesAvailable);end
    end
    fprintf(xyz,' 50 sa \n'); % new v16/v64 max acceleration to 50mm/s2 to prevent sample slipping
    fprintf(xyz,' 6 setjoybspeed \n'); % new v58  fast if joystick is pressed down
    fprintf(xyz,' 1 joystick \n'); % actually, it does not matter ??
end

if (USE_SHUTTER_CTRL)
    try
        shutter = serial(COM_SHUTTER,'BaudRate',9600,'Terminator','CR'); fopen(shutter);pause(0.5);if(shutter.BytesAvailable>0);fread(shutter,shutter.BytesAvailable);end  %SC10 thorlabs on COM6
    catch
        h=warndlg('¡ATENCIÓN! hay un problema con el controlador del shutter, suele ser fallo de la alimentación del hub USB, otras veces se ha arreglado desconectando el cable (botón 3 en el hub USB)');
        fclose(shutter);
        shutter = serial(COM_SHUTTER,'BaudRate',9600,'Terminator','CR'); fopen(shutter);pause(0.5);if(shutter.BytesAvailable>0);fread(shutter,shutter.BytesAvailable);end  %SC10 thorlabs on COM6
    end
    %to be sure the shutter is closed
    response='';
    fprintf(shutter,'ens?');pause(0.3);if(shutter.BytesAvailable>0);response=fread(shutter,shutter.BytesAvailable);end % first time a CMD_ERROR is received
    fprintf(shutter,'ens?');pause(0.3);if(shutter.BytesAvailable>0);response=fread(shutter,shutter.BytesAvailable);end
    if (numel(response) >= 6 && response(6)=='1')
        fprintf(shutter,'ens');pause(0.3);if(shutter.BytesAvailable>0);response=fread(shutter,shutter.BytesAvailable);end  % close the shutter
    end
end
if (USE_ARDUINO_QSWITCH)
    if NpulsesClean > 2  %
        CMD_QSWITCH = char('0'+NpulsesClean+2);  % parece que hay un tiempo largo antes de empezar realmente a medir, sumo 2
    else
        CMD_QSWITCH = '1';
    end
    try
        qswitch=serial(COM_QSWITCH); qswitch.BaudRate=38400;fopen(qswitch);pause(0.5);if(qswitch.BytesAvailable>0);fread(qswitch,qswitch.BytesAvailable);end  %arduino uno en com5
        %fprintf(qswitch,'0'); we would want to keep it qswitching from before
        %h=warndlg('¡ATENCIÓN! se va a encender la lámpara del láser.Asegurarse de tomar las medidas de seguridad');
        %uiwait(h);
        fprintf(qswitch,'0');
        estado_qswitch = 0;
    catch
        fclose(qswitch);
        disp('¡ATENCIÓN! ES POSIBLE QUE HAYA QUE ENCENDER EL PUERTO 1 DEL HUB USB -ARDUINO-, comprobar que está azul');
        qswitch=serial(COM_QSWITCH); qswitch.BaudRate=38400;fopen(qswitch);pause(0.5);if(qswitch.BytesAvailable>0);fread(qswitch,qswitch.BytesAvailable);end  %arduino uno en com5
        %fprintf(qswitch,'0'); we would want to keep it qswitching from before
        %h=warndlg('¡ATENCIÓN! se va a encender la lámpara del láser.Asegurarse de tomar las medidas de seguridad');
        %uiwait(h);
        fprintf(qswitch,'0');
        estado_qswitch = 0;
    end
end % if qswitch


if (USE_DUAL_PULSE_LASER)
    lotis=serial(COM_DUAL_PULSE); lotis.BaudRate=115200;fopen(lotis);pause(0.5);if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end  %arduino uno en com5
    pause(0.1); if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
    % set delays
    fwrite(lotis,['DELAY1=' num2str(DualPulse_Channel1_delay) char(13)]);
    pause(0.1); if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
    fwrite(lotis,['DELAY2=' num2str(DualPulse_Channel2_delay) char(13)]);
    pause(0.1); if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
    
    %new v50: change energy
    fwrite(lotis,['ENERGY=' num2str(round(pumpEnergy*10)) char(13)]);
    pause(0.1); if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
    h_danger=warndlg(['¡WARNING! The LOTIS-DP laser pump energy has been set to ' num2str(pumpEnergy)  'J and will be AUTOMATICALLY activated by software. Please take any safety measures from now on.']);
end % if dual_pulse

if (USE_PULSE_GENERATOR)
    pulsegen=serial(COM_PULSE_GENERATOR); pulsegen.BaudRate=38400;fopen(pulsegen);pause(0.5);if(pulsegen.BytesAvailable>0);fread(pulsegen,pulsege.BytesAvailable);end  %STM32
    pause(0.5); if(pulsegen.BytesAvailable>0);fread(pulsegen,pulsegen.BytesAvailable);end
    % set delays
    fwrite(pulsegen,[ 'q0' char(13)]); % disable qswitch pulses just in case
    estado_qswitch=0;
    pause(0.1);fwrite(pulsegen,['i' num2str(DualPulse_Channel1_delay) char(13)]);   % delay of channel 1
    pause(0.1); fwrite(pulsegen,['j' num2str(DualPulse_Channel2_delay) char(13)]);   % delay of channel 2
    pause(0.1); fwrite(pulsegen,['s' num2str(delayCaptureAvantes) char(13)]);   % delay of spectrometer
    if (UseAlternateDelay)
        pause(0.1); fwrite(pulsegen,['d' num2str(AlternateDelay) char(13)]);   % differential delay for alternate delays
        pause(0.1); fwrite(pulsegen,['a1' char(13)]);
    else
        pause(0.1); fwrite(pulsegen,['a0' char(13)]);      % no alternate delays
    end
    pause(0.1); if(pulsegen.BytesAvailable>0);fread(pulsegen,pulsegen.BytesAvailable);end
    CMD_QSWITCH = ['q' num2str(NpulsesClean+1) char(13)]; % first laser pulse
end

if (USE_LEDLIGHT)
    try
        ledlight=serial(COM_LEDLIGHT,'BaudRate',115200); fopen(ledlight);pause(0.5);if(ledlight.BytesAvailable>0);fread(ledlight,ledlight.BytesAvailable);end  %arduino en com6
        pause(0.2);
        fprintf(ledlight,'L100');
    catch
        fclose(ledlight);
        delete(instrfindall);
        
        disp('¡ATENCIÓN! no se puede acceder al puerto del controlador LED LIGHT');
        ledlight=serial(COM_LEDLIGHT,'BaudRate',115200); fopen(ledlight);pause(0.5);if(ledlight.BytesAvailable>0);fread(ledlight,ledlight.BytesAvailable);end  %arduino en com6
        pause(0.2);
        fprintf(ledlight,'L100');
    end
end % if use_ledlight


%the auto-off-exposure is not working, maybe is not enough time?
if (CAMERA)
    if (CAMERA==1)
        cameras = imaqhwinfo('avtmatlabadaptor_r2009b');
    else
        cameras = imaqhwinfo('winvideo');
    end
    if isempty(cameras.DeviceIDs)
        h=warndlg('¡ATENCIÓN! La cámara es inaccesible. Se aconseja probar: 0) Ejecutar imaqreset y después imaqtool en la linea de comandos de matlab y cerrarlo sin hacer nada, a veces resucita la cámara 1) desconectar y conectar el arduido (botón 3 en el HUB USB negro detrás de los motores. 2) salir del matlab, desconectar el cable de la cámara y volver a conectar. 3) reiniciar el ordenador.');
        return
    end
    if size(cameras.DeviceIDs,2) < 2 % webcam NOT plugged, sample camera ID =1
        ID_CAMERA_SAMPLE = 1;
    else % webcam plugged
        if (TAKE_WEBCAM_PHOTO==1)
            vid2 = videoinput('winvideo', ID_CAMERA_PEM, 'YUY2_800x600');
            src2 = getselectedsource(vid2);
            src2.WhiteBalanceMode = 'auto';
            vid2.ReturnedColorspace = 'rgb';
            triggerconfig(vid2, 'manual');
            start(vid2);
        end
        
    end
    
    if (CAMERA==1)
        try
            vid = videoinput('avtmatlabadaptor_r2009b', 1, 'F0M5_Mono8_640x480');
        catch
            %clc; clear;close all; objects = imaqfind %find video input objects in memory
            objects = imaqfind; delete(objects); %delete a video input object from memory
            vid = videoinput('avtmatlabadaptor_r2009b', 1, 'F0M5_Mono8_640x480');
        end
        src = getselectedsource(vid);
        vid.FramesPerTrigger = 1;  vid.ReturnedColorspace = 'bayer';  vid.BayerSensorAlignment = 'rggb';  src.WhiteBalanceMode = 'auto';  src.LUTMode = 'off'; % that's important
        src.GainMode = 'on';  src.FrameRate = '7';  src.ShutterMode = 'auto';  triggerconfig(vid, 'manual');
        triggerconfig(vid, 'immediate');
    else % camera 2 Amazon feb2021
        try
            vid = videoinput('winvideo', ID_CAMERA_SAMPLE, 'YUY2_1280x960');
        catch
            objects = imaqfind; delete(objects); %delete a video input object from memory
            vid = videoinput('winvideo', ID_CAMERA_SAMPLE, 'YUY2_1280x960');
        end
        src = getselectedsource(vid);
        src.WhiteBalanceMode = 'manual';
        vid.ReturnedColorspace = 'rgb';
        src.ExposureMode = 'auto'; % testing manual / auto
        src.Exposure = -2; % Amazon Camera feb2021, max value -1 changed 5nov2021 , it depends a lot on illumination 
        src.BacklightCompensation = 'off';
        triggerconfig(vid, 'manual');
        start(vid); % test v27 to speed things up
    end
    
    
    %% the laser spot affects the gain and white balance, they must be disabled
    % after a first adjust
    if (AUTOFOCUS && strcmp(AUTOFOCUS_METHOD,'laserspot')==1)
        %fwrite(theta,'l'); %% switch on the laser spot
        pause(1);
    end
    if CAMERA  frame = captureFrameGuppy;  else    frame = blankImage; end
end % if camera

if (PREVIEW)
    focusLastZ = 0; % start value for z
    lastTimeXYZ = clock; % every some seconds, the Z position is updated
    lastTimeEstTime = clock;
    scrsz=get(0,'ScreenSize');
    %h_preview=preview(vid);
    h_preview=figure('Position',[100 100 1100 900],'Name','Video Preview'); % create preview window
    hold on;
    disp('Showing video preview (flipped), click to move target, any key to close and continue...');

    esc_button_prev = uicontrol('Style', 'PushButton', 'String', 'Finish', 'Callback', 'finishPreview=1;', 'Position',[10 10 50 20], 'ForegroundColor','red');
    qswitch_button_prev = uicontrol('Style', 'PushButton', 'String', 'Activar Qsw', 'Callback', 'activar_qswitch=1;', 'Position',[70 10 70 20]);
    autofocus_button_prev = uicontrol('Style', 'PushButton', 'String', 'Autofocus', 'Callback', 'hacer_autofocus=1;', 'Position',[140 10 70 20]);
    calibraspot_button_prev = uicontrol('Style', 'PushButton', 'String', 'Calib. spot', 'Callback', 'hacer_calibracion_spotlaser=1;', 'Position',[210 10 70 20]);
    disparar_button_prev = uicontrol('Style', 'PushButton', 'String', 'Disparar', 'Callback', 'disparar=1;', 'Position',[280 10 70 20]);
    continue_button_prev = uicontrol('Style', 'PushButton', 'String', '¡Medir!', 'Callback', 'lets_go=1;', 'Position',[360 10 70 20],'ForegroundColor','green');

    left_button = uicontrol('Style', 'PushButton', 'String', '<-', 'Callback', 'fprintf(xyz,''0 0.1 0 r\n1 joystick \n'');', 'Position',[520 50 20 20]);
    right_button = uicontrol('Style', 'PushButton', 'String', '->', 'Callback', 'fprintf(xyz,''0 -0.1 0 r\n1 joystick \n'');', 'Position',[560 50 20 20]);
    superleft_button = uicontrol('Style', 'PushButton', 'String', '<-', 'Callback', 'fprintf(xyz,''0 1 0 r\n1 joystick \n'');', 'Position',[500 50 20 20]);
    superright_button = uicontrol('Style', 'PushButton', 'String', '->', 'Callback', 'fprintf(xyz,''0 -1 0 r\n1 joystick \n'');', 'Position',[580 50 20 20]);
    up_button = uicontrol('Style', 'PushButton', 'String', '^', 'Callback', 'fprintf(xyz,''-0.1 0 0 r\n1 joystick \n'');', 'Position',[540 70 20 20]);
    down_button = uicontrol('Style', 'PushButton', 'String', '|', 'Callback', 'fprintf(xyz,''+0.1 0 0 r\n1 joystick \n'');', 'Position',[540 30 20 20]);
    superup_button = uicontrol('Style', 'PushButton', 'String', '^', 'Callback', 'fprintf(xyz,''-1 0 0 r\n1 joystick \n'');', 'Position',[540 90 20 20]);
    superdown_button = uicontrol('Style', 'PushButton', 'String', '|', 'Callback', 'fprintf(xyz,''+1 0 0 r\n1 joystick \n'');', 'Position',[540 10 20 20]);
    pluz_button = uicontrol('Style', 'PushButton', 'String', '+z', 'Callback', 'fprintf(xyz,''0 0 +0.1 r\n1 joystick \n'');', 'Position',[600 70 20 20]);
    minusz_button = uicontrol('Style', 'PushButton', 'String', '-z', 'Callback', 'fprintf(xyz,''0 0 -0.1 r\n1 joystick \n'');', 'Position',[600 30 20 20]);
    superplusz_button = uicontrol('Style', 'PushButton', 'String', '++', 'Callback', 'fprintf(xyz,''0 0 1 r\n1 joystick \n'');', 'Position',[600 90 20 20]);
    superminusz_button = uicontrol('Style', 'PushButton', 'String', '--', 'Callback', 'fprintf(xyz,''0 0 -1 r\n1 joystick \n'');', 'Position',[600 10 20 20]);
        
    lessB_button = uicontrol('Style', 'PushButton', 'String', '-Brightness', 'Callback', 'stop(vid);src.Exposure = src.Exposure-1;start(vid);', 'Position',[360 50 100 20]);
    moreB_button = uicontrol('Style', 'PushButton', 'String', '+Brightness', 'Callback', 'stop(vid);if src.Exposure<-1 src.Exposure = src.Exposure+1;end;start(vid);', 'Position',[360 70 100 20]);
  
    absoluteZ_box = uicontrol('style','Edit', 'string','0', 'Position',[10 40 40 20],'backgroundcolor','w','Tag','EditField');
    absoluteZ_button = uicontrol('Style', 'PushButton', 'String', 'goto z', 'Callback', 'set_absoluteZ=1;', 'Position',[50 40 40 20]);
    
    focus_measure_text = uicontrol('style','text', 'string','Focus:', 'Position',[10 60 40 20],'Tag','EditField');
    focus_measure_value = uicontrol('style','text', 'string','0', 'Position',[50 60 50 20],'Tag','EditField');
    
    checkSurface_brightness_text = uicontrol('style','text', 'string','Brightness:', 'Position',[90 60 80 20],'Tag','EditField');
    checkSurface_brightness_value = uicontrol('style','text', 'string','??', 'Position',[160 60 30 20],'Tag','EditField');
    
    %new v58
    goto_box = uicontrol('style','Edit', 'string','0,0', 'Position',[100 40 60 20],'backgroundcolor','w','Tag','EditField');
    goto_button = uicontrol('Style', 'PushButton', 'String', 'goto x,y', 'Callback', 'pending_goto=1;', 'Position',[160 40 50 20]);
    reset_button = uicontrol('Style', 'PushButton', 'String', 'reset x,y,z=0', 'Callback', 'fprintf(xyz,'' 0 0 0 setpos \n'');', 'Position',[220 40 70 20]);
    %new v60
    xyinfo_text = uicontrol('style','text', 'string','(x,y,z) ','Position',[200 60 150 20],'Tag','EditField');

    %new v58
    stitch_x_text= uicontrol('style','text', 'string','Stitch Xmin,Xmax','Position',[710 30 100 20],'Tag','EditField');
    stitch_x_box = uicontrol('style','Edit', 'string','0,0', 'Position',[820 32 60 20],'backgroundcolor','w','Tag','EditField');
    stitch_y_text= uicontrol('style','text', 'string','Stitch Ymin,Ymax','Position',[710 10 100 20],'Tag','EditField');
    stitch_y_box = uicontrol('style','Edit', 'string','0,0', 'Position',[820 12 60 20],'backgroundcolor','w','Tag','EditField');
    set(stitch_x_box,'String',[ num2str(MinXmmStitch) ',' num2str(MaxXmmStitch)]);
    set(stitch_y_box,'String',[ num2str(MinYmmStitch) ',' num2str(MaxYmmStitch)]);
    new_stitching_button_prev = uicontrol('Style', 'PushButton', 'String', 'New stitch', 'Callback', 'do_new_stitching=1;', 'Position',[650 80 70 20]);
    load_stitching_button_prev = uicontrol('Style', 'PushButton', 'String', 'Load stitch', 'Callback', 'do_load_stitching=1;', 'Position',[650 60 70 20]);
    new_polygon_button_prev = uicontrol('Style', 'PushButton', 'String', 'New polygon', 'Callback', 'do_new_polygon=1;', 'Position',[730 80 70 20]);
    load_polygon_button_prev = uicontrol('Style', 'PushButton', 'String', 'Load polygon', 'Callback', 'do_load_polygon=1;', 'Position',[730 60 70 20]);
    new_linepath_button_prev = uicontrol('Style', 'PushButton', 'String', 'New linePath', 'Callback', 'do_new_linepath=1;', 'Position',[810 80 70 20]);
    load_linepath_button_prev = uicontrol('Style', 'PushButton', 'String', 'Load linePath', 'Callback', 'do_load_linepath=1;', 'Position',[810 60 70 20]);

    %new v13/v64
    description_text= uicontrol('style','text', 'string','Description:','Position',[220 190 90 20],'Tag','EditField','HorizontalAlignment','right');
    description_box = uicontrol('style','Edit', 'string',EXPERIMENT_DESCRIPTION, 'Position',[320 192 250 20],'backgroundcolor','w','Tag','EditField','HorizontalAlignment','left');
    path_text= uicontrol('style','text', 'string','Path:','Position',[220 170 90 20],'Tag','EditField','HorizontalAlignment','right');
    path_box = uicontrol('style','Edit', 'string',EXPERIMENT_NAME, 'Position',[320 172 200 20],'backgroundcolor','w','Tag','EditField','HorizontalAlignment','left');
    experiment_text= uicontrol('style','text', 'string','Subpath:','Position',[220 150 90 20],'Tag','EditField','HorizontalAlignment','right');
    experiment_box = uicontrol('style','Edit', 'string',SUBEXPERIMENT_NAME, 'Position',[320 152 250 20],'backgroundcolor','w','Tag','EditField','HorizontalAlignment','left');
    npulses_text= uicontrol('style','text', 'string','Npulses:','Position',[720 190 90 20],'Tag','EditField','HorizontalAlignment','right');
    npulses_box = uicontrol('style','Edit', 'string',Npulses, 'Position',[820 192 60 20],'backgroundcolor','w','Tag','EditField','HorizontalAlignment','left');
    step_text= uicontrol('style','text', 'string','Step(mm):','Position',[720 170 90 20],'Tag','EditField','HorizontalAlignment','right');
    step_box = uicontrol('style','Edit', 'string',param1Step, 'Position',[820 172 80 20],'backgroundcolor','w','Tag','EditField','HorizontalAlignment','left');
    estimated_time_text = uicontrol('style','text', 'string','Estimated time: ?','Position',[750 130 200 20],'Tag','EditField','HorizontalAlignment','left');
    type1D2D_text= uicontrol('style','text', 'string','Type: 1D (deltaR)','Position',[270 130 200 20],'Tag','EditField','HorizontalAlignment','left');
    if strcmp(param1Name,'deltaR') 
         type1D2D_text.String = 'Type: 1D (deltaR)';
         absoluteMovement = 0;
    elseif strcmp(param1Name,'posX') || strcmp(param1Name,'posY')
         type1D2D_text.String = 'Type: 2D (posX,posY)';
         absoluteMovement = 1;
    else
         type1D2D_text.String = ['Type: ' param1Name];
    end

    finishPreview=0;
    lets_go=0;
    disparar=0;
    hacer_autofocus=0;
    hacer_calibracion_spotlaser=0;
    do_new_stitching=0;
    do_load_stitching=0;
    do_new_polygon=0;
    do_load_polygon=0;
    do_new_linepath=0;
    do_load_linepath=0;
    set_absoluteZ=0; % changed to 1 by the set absZ button
    k=[]; set(h_preview,'keypress','k=get(gcf,''currentchar'');');
    previewTime = tic;
    frame = captureFrameGuppy;
    h_ims_prev=imshow(frame,'InitialMagnification', 67);
    %v64 target() drawn only once
    target(pixel_shotX,pixel_shotY); % draw a cross around
    if (ESCALA==1) % draw 1mm ticks starting at the laser spot
        if calibrated_pixels_mm_XY_X > 300   % new v58
            mm1o01=calibrated_pixels_mm_XY_X*0.1; % line every 100um
            numberOflines=7;
        else
            mm1o01=calibrated_pixels_mm_XY_X; % line every 1mm
            numberOflines=3;
        end
        %hold on;
        for horticks=1:numberOflines % cortos hacia la izquierda
            line([pixel_shotX-(horticks-1)*mm1o01  pixel_shotX-(horticks-1)*mm1o01],[0 300],'Color','w');
            line([pixel_shotX-(horticks-1)*mm1o01+1  pixel_shotX-(horticks-1)*mm1o01+1],[0 300],'Color','k');
        end
        for horticks=1:numberOflines  % largos hacia la derecha
            line([pixel_shotX+(horticks)*mm1o01  pixel_shotX+(horticks)*mm1o01],[0 600],'Color','w');
            line([pixel_shotX+(horticks)*mm1o01+1  pixel_shotX+(horticks)*mm1o01+1],[0 600],'Color','k');
        end
    end
    pending_goto=0;
    while 1
        %feature('memstats');
        frame = captureFrameGuppy;
        if (~isempty(k) || ~ishandle(h_preview) || finishPreview == 1 || lets_go==1)  % if it has been closed
            break;
        end
        
        %new v55 brightness value around laser shot (mean of 5x5 pixels)
        matrixSpot = frame(pixel_shotY-2:pixel_shotY+2,pixel_shotX-2:pixel_shotX+2);
        brightnessSpot = mean(mean(matrixSpot));
        set(checkSurface_brightness_value,'String',num2str(brightnessSpot,'%3.0f'));
        set(h_ims_prev,'CData',frame); % new! v28 update image without a new imshow, promise speed things up
        if (0)
            if (ESCALA==1) % draw 1mm ticks starting at the laser spot
                if calibrated_pixels_mm_XY_X > 300   % new v58
                    mm1o01=calibrated_pixels_mm_XY_X*0.1; % line every 100um
                    numberOflines=7;
                else
                    mm1o01=calibrated_pixels_mm_XY_X; % line every 1mm
                    numberOflines=3;
                end
                hold on;
                for horticks=1:numberOflines % cortos hacia la izquierda
                    line([pixel_shotX-(horticks-1)*mm1o01  pixel_shotX-(horticks-1)*mm1o01],[0 300],'Color','w');
                    line([pixel_shotX-(horticks-1)*mm1o01+1  pixel_shotX-(horticks-1)*mm1o01+1],[0 300],'Color','k');
                end
                for horticks=1:numberOflines  % largos hacia la derecha
                    line([pixel_shotX+(horticks)*mm1o01  pixel_shotX+(horticks)*mm1o01],[0 600],'Color','w');
                    line([pixel_shotX+(horticks)*mm1o01+1  pixel_shotX+(horticks)*mm1o01+1],[0 600],'Color','k');
                end
            end % if ESCALA
            target(pixel_shotX,pixel_shotY); % draw a cross around
        end % do now draw target nor lines
        if (activar_qswitch == 1 && (USE_ARDUINO_QSWITCH || USE_DUAL_PULSE_LASER || USE_PULSE_GENERATOR)) % activates arduino relay or lotis shutter&qswitch
            activar_qswitch = 0;
            if (estado_qswitch==0)
                estado_qswitch=1;
                if (USE_ARDUINO_QSWITCH==1)
                    fprintf(qswitch,'1'); % activates relay
                end
                if (USE_DUAL_PULSE_LASER)
                    fwrite(lotis,['SHUT1ON' char(13)]);
                    fwrite(lotis,['QSW1ON' char(13)]);
                    if (DualPulse) % also channel 2
                        fwrite(lotis,['SHUT2ON' char(13)]);
                        fwrite(lotis,['QSW2ON' char(13)]);
                    end
                    pause(0.1);
                    if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
                end
                if (USE_PULSE_GENERATOR==1)
                    fprintf(pulsegen,['q1' char(13)]); % activates qswitch pulses with no 'noise' pulses
                end
                set(qswitch_button_prev , 'string','Desactivar Qsw');
                drawnow;
            else
                estado_qswitch=0;
                if (USE_ARDUINO_QSWITCH==1)
                    fprintf(qswitch,'0');
                end
                if (USE_DUAL_PULSE_LASER)
                    fwrite(lotis,['QSW1OFF' char(13)]);
                    fwrite(lotis,['SHUT1OFF' char(13)]);
                    if (DualPulse) % also channel 2
                        fwrite(lotis,['QSW2OFF' char(13)]);
                        fwrite(lotis,['SHUT2OFF' char(13)]);
                    end
                    pause(0.1);
                    if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
                end
                if (USE_PULSE_GENERATOR==1)
                    fprintf(pulsegen,['q0' char(13)]); % remove qswitch pulses
                end
                set(qswitch_button_prev , 'string','Activar Qsw');
                drawnow;
            end
        end        
        
        if MOTORS && etime(clock,lastTimeXYZ) > 2 % every 3 seconds the Z position is updated,  and the focus index
            fprintf(xyz,' p \r1 j'); % get position and enable joystick again
            pause(0.05);
            current_coords = fscanf(xyz,'%f %f %f'); % this reads the current coordinates x y z
            absX = current_coords(2);
            absY = current_coords(1);
            absZ = current_coords(3);
            %set(absoluteZ_box,'String',absZ); % z value is already printed elsewhere
            infoxy=['(x,y,z)=(',num2str(absX,'%.2f') ,',',num2str(absY,'%.2f'),',',num2str(absZ,'%.2f'),')'];
            set(xyinfo_text,'String',infoxy);
            frame = captureFrameGuppy; % just in case there is a refresh issue
            bw = rgb2gray(frame);
            f =  fmeasure(bw,'GRAS',[pixel_shotX-100 pixel_shotY-100 100 100]);
            set(focus_measure_value,'String',num2str(f,'%.2f'));
            lastTimeXYZ = clock;
        end
        
        if MOTORS && set_absoluteZ==1 % change Z position new! v52
            fprintf(xyz,' p \r');
            pause(0.05);
            current_coords = fscanf(xyz,'%f %f %f'); % this reads the current coordinates x y z
            absX = current_coords(2);
            absY = current_coords(1);
            absZ = current_coords(3);
            absZnow = str2double(absoluteZ_box.String); % user value
            moveString = sprintf(' %0.3f %0.3f %0.3f m \r',absY, absX, absZnow); % absolute values, command move (m)
            fprintf(xyz,moveString); %move axis
            set_absoluteZ=0;
        end
        if MOTORS && hacer_autofocus==1 % performs an initial autofocus
            set(autofocus_button_prev,'Enable','off')
            disp('Doing autofocus...');     % el archivo experiment.txt no está abierto
            % testing v52
            %i=1;
            %for n=-2:0.1:2
            %    val = minimizaFocus(n);
            %    disp([n val]);
            %    focos_z(i)=n;
            %    focos_f(i)=val;
            %    i=i+1;
            %end
            %apuntamos la focalización antes de empezar
            focus_antes = minimizaFocus(0);
            options = optimset('Display','notify','TolX',0.001,'MaxFunEvals',20);
            [z focusVal exitFlag] = fminbnd(@minimizaFocus,-3,3,options); % minimiza en un rango RELATIVO de +/-2 milímetros, el rango es pequeño porque con rangos de +/-5mm se quedaba en mínimos locales ¿?¿
            if (exitFlag == 0) % no se ha encontrado foco
                disp('Focus not found!');
                sound(sound0);
            end
            minimizaFocus(z); % move to the best focus
            focusString = ['Focus index before=' num2str(-focus_antes) ' after=' num2str(-focusVal) ' relative defocus z=' num2str(z) ' absolute position z=(see box in window)' ];
            disp(focusString);     % el archivo experiment.txt no está abierto
            set(autofocus_button_prev,'Enable','on')
            hacer_autofocus=0;
            fprintf(xyz,' 1 joystick \n'); % enable joystick
        end
        %new v33 change the calibration of the laser spot position
        if(hacer_calibracion_spotlaser==1)
            %frame_antes = captureFrameGuppy;
            h_cal=warndlg('Calibración: disparar alguna vez con modo INTERNO del laser y pulsar el botón para posicionar la cruceta');
            uiwait(h_cal);
            %frame_despues =  captureFrameGuppy;
            %h_click=figure('Name','Laser spot position calibration - click in the center of the laser spot'); % create new large window
            %h_ims_prev=imshow(frame,'InitialMagnification', 67);
            %frame_diff = frame_despues-frame_antes;
            %frame_diff = frame_diff ./ max(max(max(frame_diff)));
            %imshow(frame_despues-frame_antes); hold on;
            frame = captureFrameGuppy;
            set(h_ims_prev,'CData',frame); % new! v28 update image without a new imshow, promise speed things up
            [xc, yc] = ginput(1);
            xc=round(xc);yc=round(yc);
            % Construct a questdlg with three options
            sino= questdlg(['Last: (' num2str(pixel_shotX,'%4.0f') ',' num2str(pixel_shotY,'%4.0f') ')  New: (' num2str(xc,'%4.0f') ',' num2str(yc,'%4.0f') ') Change?' ], 'Yes', 'No');
            switch sino
                case 'Yes'
                    pixel_shotX = xc; pixel_shotY = yc;
                    save('last_laserspot_coordinates.mat','pixel_shotX','pixel_shotY')
                case 'No'
                    % do nothing
            end
            hacer_calibracion_spotlaser=0;
            %close(h_click);
        end
        %new v58 if pending_goto=1
        if MOTORS && pending_goto==1
            fprintf(xyz,' p \r');
            pause(0.05);
            current_coords = fscanf(xyz,'%f %f %f'); % this reads the current coordinates x y z
            goto_coords = sscanf(goto_box.String,'%f,%f'); % this reads the current coordinates x y z
            % move X,Y (positions 2 ans 1), but not in Z
            moveString = sprintf(' %0.3f %0.3f %0.3f m \r',goto_coords(2), goto_coords(1), current_coords(3)); % absolute values, command move (m)
            fprintf(xyz,moveString); %move axis
            pause(0.05);
            fprintf(xyz,' 1 joystick \n'); % enable joystick
            pending_goto=0;
        end
        
        %new v63 load image
        if do_load_stitching
            [ifile,ipath]= uigetfile('*.fig','Load figure with stitched image',filenameImageStitched);
            if ~isequal(ifile,0)
               filenameImageStitched = strcat(ipath,ifile);
               try
                   h_stitch=openfig(filenameImageStitched);
                   figure(h_stitch);
                   disp(['Loaded ' filenameImageStitched]);
                   StitchedIsLoadedNotCreated=1;
                   haveImageStitched=1;
                   p=findobj(h_stitch,'Type','images.roi.polygon'); % check if polygon is already drawn
                   if ~isempty(p) % it includes a polygon, so let's load it and change some things...
                       roiPol = drawpolygon('Position',p(1).Position);
                       delete(p); % new v64
                       set(roiPol,'Color','green');
                       havePolygon=1;
                       param1Name='posX';
                       param2Name='posY';
                       MANUAL_TRACKING=0;
                       type1D2D_text.String = 'Type: 2D (posX,posY)';
                       absoluteMovement = 1;
                   end
                   p=findobj(h_stitch,'Type','images.roi.polyline'); % check if linepath is already drawn
                   if ~isempty(p) % there is a linepath, if both linepath and polygon are in the image, this takes precedence
                       roiLine = drawpolyline('Position',p(1).Position);
                       delete(p); % new v64
                       set(roiLine,'Color','green');
                       haveLinePath=1;
                       param1Name='deltaR';
                       param2Name='none';
                       param2From = 0;
                       param2To = 0;
                       MANUAL_TRACKING=1;
                       type1D2D_text.String = 'Type: 1D (deltaR)';
                       absoluteMovement = 0;
                   end
               catch
                   error('ERROR!! the filename with the stitchedImage cannot be read, check name');
               end
            end
            do_load_stitching=0;
        end
 
        %new v63 new polygon
        if do_new_polygon && exist('h_stitch','var')
            if ishandle(h_stitch) % still, could have a closed figure
                figure(h_stitch);
                roiPol = drawpolygon; % https://es.mathworks.com/help/images/ref/drawpolygon.html
                havePolygon=1;
                do_new_polygon=0;
                param1Name='posX';
                param2Name='posY';
                MANUAL_TRACKING=0;
                type1D2D_text.String = 'Type: 2D (posX,posY)';
                absoluteMovement = 1;
            end
        end
        %new v63 load polygon
        if do_load_polygon && exist('h_stitch','var')
            [pfile,ppath]= uigetfile('*.mat','Load polygon',filenamePolygon);
            if ~isequal(pfile,0)
                filenamePolygon= strcat(ppath,pfile);
                try
                    load(filenamePolygon);
                catch
                    disp('ERROR!! the filename with the Polygon cannot be read, check name');
                end
                try
                    disp(['Loaded ' filenamePolygon]);
                    havePolygon=1;
                    figure(h_stitch);
                    roiPol = drawpolygon('Position',roiPol.Position); % print again the loaded polygon
                    param1Name='posX';
                    param2Name='posY';
                    MANUAL_TRACKING=0;
                    type1D2D_text.String = 'Type: 2D (posX,posY)';
                    absoluteMovement = 1;
                catch
                    disp('Warning!! there is no stitchedImage, but polygon was loaded');
                end
            end
            do_load_polygon=0;
        end
        %new v63 new linepath
        if do_new_linepath && exist('h_stitch','var')
            if ishandle(h_stitch) % still, could have a closed figure
                figure(h_stitch);
                roiLine = drawpolyline; % https://es.mathworks.com/help/images/ref/drawpolygon.html
                haveLinePath=1;
                do_new_linepath=0;
                param1Name='deltaR';
                param2Name='none';
                param2From = 0;
                param2To = 0;
                MANUAL_TRACKING=1;
                type1D2D_text.String = 'Type: 1D (deltaR)';
                absoluteMovement = 0;
            end
        end
        %new v63 load linepath
        if do_load_linepath && exist('h_stitch','var')
            [lfile,lpath]= uigetfile('*.mat','Load linePath',filenameLinePath);
            if ~isequal(lfile,0)
                filenameLinePath= strcat(lpath,lfile);
                try
                    load(filenameLinePath);
                catch
                    disp('ERROR!! the filename with the linePath cannot be read, check name');
                end
                try
                    disp(['Loaded ' filenameLinePath]);
                    haveLinePath=1;
                    figure(h_stitch);
                    roiLine = drawpolyline('Position',roiLine.Position); % print again the loaded linePath
                    param1Name='deltaR';
                    param2Name='none';
                    param2From = 0;
                    param2To = 0;
                    MANUAL_TRACKING=1;
                    type1D2D_text.String = 'Type: 1D (deltaR)';    
                    absoluteMovement = 0;
                catch
                    disp('Warning!! there is no stitchedImage, but linePath was loaded');
                end
            end
            do_load_linepath=0;
        end



        %new v60 stich an image of larger area
        if do_new_stitching==1
            do_new_stitching=0;
            fprintf('Creating stitched image...');
            ss = sscanf(stitch_x_box.String,'%f,%f'); % get stitching corners from text boxes
            MinXmmStitch = ss(1);MaxXmmStitch = ss(2);
            ss = sscanf(stitch_y_box.String,'%f,%f');
            MinYmmStitch = ss(1);MaxYmmStitch = ss(2);
            if MinXmmStitch>MaxXmmStitch
                tmpvar=MinXmmStitch;MinXmmStitch=MaxXmmStitch;MaxXmmStitch=tmpvar;
            end
            if MinYmmStitch>MaxYmmStitch
                tmpvar=MinYmmStitch;MinYmmStitch=MaxYmmStitch;MaxYmmStitch=tmpvar;
            end
            % vinCorr=im2double(imread('vignetting_correction.png')); % trying to correct vignetting REMOVED v64 it does not work well
            fprintf(xyz,' p \r');    pause(0.05);
            starting_coords = fscanf(xyz,'%f %f %f'); % this reads the current coordinates x y z
            stWidth = ceil(dxStitch*calibrated_pixels_mm_XY_X);
            stHeight = ceil(dyStitch*calibrated_pixels_mm_XY_Y);
            stPixelsX = size(MinXmmStitch:dxStitch:MaxXmmStitch,2)*stWidth;
            stPixelsY = size(MinYmmStitch:dyStitch:MaxYmmStitch,2)*stHeight;
            stitchImage  = zeros(stPixelsY,stPixelsX,3,'uint8');
            h_stitch = figure('Name','Image stitching'); % 'Position',[100 100 1000 800]
            finishStitch=0;
            esc_button_stitch = uicontrol('Style', 'PushButton', 'String', 'Finish', 'Callback', 'finishStitch=1;', 'Position',[10 10 50 20], 'ForegroundColor','red');
            %capture_polygon_stitch = uicontrol('Style', 'PushButton', 'String', 'Capture Polygon', 'Callback', 'do_capture_polygon=1;', 'Position',[70 10 90 20]);
            %capture_linePath_stitch = uicontrol('Style', 'PushButton', 'String', 'Capture Measuring Path', 'Callback', 'do_capture_linePath=1;', 'Position',[170 10 120 20]);
            %remove_artifacts_stitch = uicontrol('Style', 'PushButton', 'String', 'Remove artifacts', 'Callback', 'do_remove_artifacts=1;', 'Position',[300 10 90 20]);
            do_capture_polygon=0;
            do_capture_linePath=0;
            do_remove_artifacts=0;
            
            %hsi = image([ MinXmmStitch-dxStitch/2  MaxXmmStitch+dxStitch/2] , [MinYmmStitch-dyStitch/2 MaxYmmStitch+dyStitch/2],flipud(stitchImage));
            hsi = imshow(stitchImage);
            pbaspect([abs(MaxXmmStitch-MinXmmStitch) abs(MaxYmmStitch-MinYmmStitch) 1]) % this is the real aspect ratio
            xOffsetToCenter = pixel_shotX - 1280/2; % pixels from laser spot to image center
            yOffsetToCenter = pixel_shotY - 960/2; 
            for stloopy=MinYmmStitch:dyStitch:MaxYmmStitch
                for stloopx=MinXmmStitch:dxStitch:MaxXmmStitch
                    fprintf(xyz,' p \r'); pause(0.05);
                    current_coords = fscanf(xyz,'%f %f %f'); % this reads the current coordinates x y z
                    distance=sqrt((current_coords(1)-stloopy)^2+(current_coords(2)-stloopx)^2);
                    moveString = sprintf(' %0.3f %0.3f %0.3f m \r',stloopy-yOffsetToCenter/calibrated_pixels_mm_XY_Y, stloopx+xOffsetToCenter/calibrated_pixels_mm_XY_X, starting_coords(3)); % absolute values, command move (m)
                    fprintf(xyz,moveString); %move axis
                    %new v16/v64 wait for motor stop
                    beforeMovingForStitching=tic();
                    statusXYZ=1;
                    while statusXYZ==1        %is moving,
                        pause(0.05);
                        if(xyz.BytesAvailable>0);fread(xyz,xyz.BytesAvailable);end
                        fprintf(xyz,' status \r');
                        statusXYZ = fscanf(xyz,'%f'); % status word
                    end
                    %disp([' time to motor stop in stitching: ' num2str(toc(beforeMovingForStitching)) ' distance: ' num2str(distance) ]);
                    pause(0.6); % camera is set to 7fps, this should get a completely-stopped image, I hope , changed to 0.8 17jul2024 because edge mistmatch were too visible, back to 0.6, it does not make any differnce
                    coordx = stloopx + xOffsetToCenter;
                    coordy = stloopy + yOffsetToCenter;
                    frame=captureFrameGuppy;                    
                    % vignetting correction using a defocused image REMOVED v16,it does not work well in last tests
                    %fraVigCorr = im2double(frame)./vinCorr;
                    %fraVigCorr(fraVigCorr>1.0)=1.0; fraVigCorr = uint8(fraVigCorr*255);% clip at 1 and back to uint8 
                    %fraVigCorr = rescale(fraVigCorr, 0, 255,'InputMax',1);fraVigCorr = uint8(fraVigCorr); % not supported in matlab 2016
                    %stRectangle = fraVigCorr(960/2-stHeight/2:960/2+stHeight/2,1280/2-stWidth/2:1280/2+stWidth/2,:); %get rectangle around target with vignetting correction
                    stRectangle = frame(round(960/2-stHeight/2):round(960/2+stHeight/2),round(1280/2-stWidth/2):round(1280/2+stWidth/2),:); %get rectangle around target WITHOUT vignetting correction
                    
                    stPixelOffsetX = round((stloopx-MinXmmStitch)*calibrated_pixels_mm_XY_X+1);
                    stPixelOffsetY = round(stPixelsY - stHeight - (stloopy-MinYmmStitch)*calibrated_pixels_mm_XY_Y+1); % image coordinates start upper-left corner
                    stitchImage(stPixelOffsetY:stPixelOffsetY+stHeight,stPixelOffsetX:stPixelOffsetX+stWidth,:) = stRectangle; % copy rectangle to global image
                    % add cross and coordinates
                    %target(stPixelOffsetX+stWidth/2,stPixelOffsetY+stHeight/2); % great cross but can losses focus
                    % UNCOMMENT these two lines to draw the cross and print the coordinates
                    %stitchImage = insertShape(stitchImage,'line', [ stPixelOffsetX+round(stWidth/2)-10 stPixelOffsetY+round(stHeight/2) stPixelOffsetX+round(stWidth/2)+10 stPixelOffsetY+round(stHeight/2)  ;  stPixelOffsetX+round(stWidth/2) stPixelOffsetY+round(stHeight/2)-10 stPixelOffsetX+round(stWidth/2) stPixelOffsetY+round(stHeight/2)+10  ], Color='black', LineWidth=3); % couldn't make it work
                    %stitchImage = insertText(stitchImage, [ stPixelOffsetX+round(stWidth/2)+5 stPixelOffsetY+round(stHeight/2)+5 ], [ num2str(stloopx) ',' num2str(stloopy) ], FontSize=24, BoxOpacity=0, TextColor='black' );
                    if ( ~ishandle(h_stitch) || finishStitch == 1 )  % if it has been closed or wanted to leave
                        break;
                    end
                    set(hsi,'CData',stitchImage);
                    %figure(h_stitch);
                    %hsi = imshow(stitchImage,'InitialMagnification', 50);
                end
                if ( ~ishandle(h_stitch) || finishStitch == 1 )  % if it has been closed or wanted to leave
                    break;
                end
            end
            % clear vinCorr; % new v11 takes a lot of memory REMOVED v16
            if ishandle(h_stitch)
                haveImageStitched=1;
                disp('done!');
                figure(h_stitch);
                hsi = image([ MinXmmStitch-dxStitch/2  MinXmmStitch+size(MinXmmStitch:dxStitch:MaxXmmStitch,2)*dxStitch-dxStitch/2] , [MinYmmStitch-dyStitch/2 MinYmmStitch+size(MinYmmStitch:dyStitch:MaxYmmStitch,2)*dyStitch-dyStitch/2],flipud(stitchImage));
                set(gca,'YDir','normal');
                pbaspect([abs(MaxXmmStitch-MinXmmStitch) abs(MaxYmmStitch-MinYmmStitch) 1]) % this is the real aspect ratio
            else
                disp('cancelled!');
            end
            moveString = sprintf(' %0.3f %0.3f %0.3f m \r',starting_coords(1), starting_coords(2), starting_coords(3)); % absolute values, command move (m)
            fprintf(xyz,moveString); %move back to where it was
            % new v60 interactive tools with stitched image
            while(0)
                drawnow;
                if ~ishandle(h_stitch)  % if it has been closed or wanted to leave
                    haveImageStitched=0; % the window is closed
                    break;
                end
                if finishStitch == 1  %just leave to preview
                    break;
                end
                if do_capture_polygon % capture a polygonal ROI for 2d measurements
                    roiPol = drawpolygon; % https://es.mathworks.com/help/images/ref/drawpolygon.html
                    havePolygon=1;
                    do_capture_polygon=0; % but keep in the loop, just in case you want another try
                end
                if do_capture_linePath
                    roiLine = drawpolyline;
                    haveLinePath=1;
                    do_capture_linePath=0;
                end
                if do_remove_artifacts
                    do_remove_artifacts=0;
                    fprintf('Removing artifacts, could take a while and a LOT of memory and CPU ... ');
                    stitchImageOriginal = stitchImage;
                    Y = im2double( stitchImage );
                    lambda = 1.0; tol = 0; nite = 20;
                    [L, E] = imdecpattern( Y, lambda, tol, nite ); %https://es.mathworks.com/matlabcentral/fileexchange/63145-pattern-artifact-removal-from-an-image?s_tid=blogs_rc_6
                    stitchImage=im2uint8(L); % the 'decomposed' image
                    set(hsi,'CData',flipud(stitchImage));
                    %set(gca,'YDir','normal');
                    disp('done!');
                end
            end
            figure(h_preview);
        end
        
        %new v13/v64 calcs estimated measuring time for 2D polygons
        try
            if havePolygon && etime(clock,lastTimeEstTime) > 2 % every 2 seconds, only if polygon is defined
                xmin=min(roiPol.Position(:,1));
                xmax=max(roiPol.Position(:,1));
                ymin=min(roiPol.Position(:,2));
                ymax=max(roiPol.Position(:,2));
                param1Step=str2double(step_box.String);
                param2Step=param1Step;
                ejex = xmin:param1Step:xmax;
                ejey = ymin:param2Step:ymax;
                [gridx,gridy]=meshgrid(ejex,ejey);
                effectivePoints =sum(inROI(roiPol,gridx,gridy));
                timePerPoint = (str2double(npulses_box.String) *0.1) + 2.33 ; % the npulsesx0.1secs plus "overhead"
                estimated_time_text.String = ['Estimated time: ' datestr(seconds(timePerPoint*effectivePoints),'DD:HH:MM:SS') ' (' num2str(Npixeles*Npulses*Nparam1*Nparam2*4/1024/1024/1024,'%2.0f') 'GB)'];
                lastTimeEstTome = clock;
            end
        catch
            estimated_time_text.String = ['Estimated time: unkown']; %polygon can be deleted
        end    
        
            %new v38 display the avantes spectrum in real time
            %new v54 with absolute energy and line heights. subplots: 1-3 avantes; 4-6 pimax; 7 energy bar :)
            if (disparar==1)
                finishDisparar=0;
                disparar_save_spectrum=0;
                disparar_dualpulse=0;
                disparar_delaycapture=0;
                
                % display spectrum
                energia = 10; %relative energy of spectrum
                energia_abs = 0; % absolute (sum of raw values)
                h_disparar = figure('Position',[ 100 100 scrsz(3)/2 scrsz(4)/1.3 ] ,'Name','laser shots to align the fiber');
                
                hold off;
                %axes('position',[.1 .1 .8 .6]);
                subplot(9,1,1:3,'Parent',h_disparar);
                h_plot_disparar = plot(avantes_lambdas,avantes_lambdas);
                axis([ 150 950 0 50000]);
                %axes('position',[.1 .7 .8 .2]);
                if (USE_REMOTE_PIMAX)
                    subplot(9,1,4:6,'Parent',h_disparar);
                    %h_plot_rempimax = plot(avantes_lambdas,avantes_lambdas); %not
                    %known yet
                    %axis([ 150 950 0 66000]);
                    sock = msconnect(remote_pimax_ip,3000);
                    comando = ['npulses=1'];
                    mssend(sock, comando, 2);ack = msrecv(sock,3); msclose(sock);
                end
                subplot(9,1,7,'Parent',h_disparar);
                h_plot_disparar_barra = barh(energia);
                xlim([0 100]);
                
                subplot(9,1,8:9,'Parent',h_disparar);
                h_plot_disparar_lines = barh(lineMonHeights);
                set(gca,'XScale','log');
                xlim([800 65000]); % new v58 raw value to better detect saturation
                yticklabels(lineMonitoring_name);
                
                %subplot(9,1,9); % bottom part - interestingly, it creates a new plot
                disparar_esc_button = uicontrol('Style', 'PushButton', 'String', 'Finish', 'Callback', 'finishDisparar=1;', 'Position',[100 10 50 20],'ForegroundColor','red');
                disparar_save_spectrum_button = uicontrol('Style', 'PushButton', 'String', 'Save spectrum', 'Callback', 'disparar_save_spectrum=1;', 'Position',[10 10 80 20]);
                disparar_dualpulse_button = uicontrol('Style', 'PushButton', 'String', 'set DP delay', 'Callback', 'disparar_164000=1;', 'Position',[300 10 80 20]);
                disparar_text_dualpulse = uicontrol('style','Edit', 'string','5000',  'Position',[200 10 80 20], 'backgroundcolor','w', 'Tag','EditField');
                disparar_delaycapture_button = uicontrol('Style', 'PushButton', 'String', 'set spec delay', 'Callback', 'disparar_var=1;', 'Position',[550 10 80 20]);
                disparar_text_delaycapture = uicontrol('style','Edit', 'string','1000', 'Position',[450 10 80 20],'backgroundcolor','w','Tag','EditField');
                disparar_lines = uicontrol('style','text', 'string','line heights=?', 'Position',[10 45 1100 20]);
                %new v60 not yet implemented
                %disparar_xyinfo_text = uicontrol('style','text', 'string','(X,Y) ','Position',[280 60 90 20],'Tag','EditField');
                %new v61
                left_button = uicontrol('Style', 'PushButton', 'String', '<-', 'Callback', 'fprintf(xyz,''0 0.1 0 r\n1 joystick \n'');', 'Position',[720 140 20 20]);
                right_button = uicontrol('Style', 'PushButton', 'String', '->', 'Callback', 'fprintf(xyz,''0 -0.1 0 r\n1 joystick \n'');', 'Position',[760 140 20 20]);
                up_button = uicontrol('Style', 'PushButton', 'String', '^', 'Callback', 'fprintf(xyz,''-0.1 0 0 r\n1 joystick \n'');', 'Position',[740 160 20 20]);
                down_button = uicontrol('Style', 'PushButton', 'String', '|', 'Callback', 'fprintf(xyz,''+0.1 0 0 r\n1 joystick \n'');', 'Position',[740 120 20 20]);
                plusz_button = uicontrol('Style', 'PushButton', 'String', '+z', 'Callback', 'fprintf(xyz,''0 0 +0.1 r\n1 joystick \n'');', 'Position',[800 160 20 20]);
                minus_button = uicontrol('Style', 'PushButton', 'String', '-z', 'Callback', 'fprintf(xyz,''0 0 -0.1 r\n1 joystick \n'');', 'Position',[800 120 20 20]);
                superplusz_button = uicontrol('Style', 'PushButton', 'String', '++', 'Callback', 'fprintf(xyz,''0 0 1 r\n1 joystick \n'');', 'Position',[780 160 20 20]);
                superminusz_button = uicontrol('Style', 'PushButton', 'String', '--', 'Callback', 'fprintf(xyz,''0 0 -1 r\n1 joystick \n'');', 'Position',[780 120 20 20]);
                
                % activate q-switch
                if (USE_ARDUINO_QSWITCH) estado_qswitch=0; fprintf(qswitch,'1');    end
                % open shutter
                if (USE_SHUTTER_CTRL)
                    fprintf(shutter,'ens');pause(0.2);if(shutter.BytesAvailable>0);response=fread(shutter,shutter.BytesAvailable);end  % open the shutter
                end
                if (USE_DUAL_PULSE_LASER)
                    fwrite(lotis,['SHUT1ON' char(13)]);
                    fwrite(lotis,['QSW1ON' char(13)]);
                    if (DualPulse) % also channel 2
                        fwrite(lotis,['SHUT2ON' char(13)]);
                        fwrite(lotis,['QSW2ON' char(13)]);
                    end
                    pause(0.1);
                    if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
                end
                if (USE_PULSE_GENERATOR==1)
                    fprintf(pulsegen,'q1'); % activates qswitch pulses with no 'noise' pulses
                end
                %do get the spectrum
                medido_remote=0;
                cuando_remote = now();
                while(true)
                    if(USE_AVANTES) % collect Npulse spectra from Avantes in hardware-external-trigger
                        spectrometer('measure',h1,1); spectrometer('measure',h2,1);  spectrometer('measure',h3,1); spectrometer('measure',h4,1);  %
                        spectrometer('measure',h5,1); spectrometer('measure',h6,1);  spectrometer('measure',h7,1); spectrometer('measure',h8,1);  %
                    end
                    if (USE_REMOTE_PIMAX && medido_remote==0)
                        sock = msconnect(remote_pimax_ip,3000);
                        mssend(sock,'measure',2); % 2=timeout
                        cuando_remote = now();
                        medido_remote = 1;
                    end
                    if (USE_AVANTES)
                        ready=false;
                        while (ready==false)
                            if ready==false
                                ready=spectrometer('ready',h8); % h5 es el master, aunque con el cable pin6 da igual
                            end
                            pause(0.001); %seconds !!
                        end
                        pause(0.05); % just in case not all channels are ready
                        myData1 = spectrometer('getdata',h1);myData2 = spectrometer('getdata',h2);myData3 = spectrometer('getdata',h3);myData4 = spectrometer('getdata',h4);myData5 = spectrometer('getdata',h5);myData6 = spectrometer('getdata',h6);myData7 = spectrometer('getdata',h7);myData8 = spectrometer('getdata',h8);
                        %all_data = cat(1,myData1(1:1921),myData2(1:1867),myData3(1:1786),myData4(1:1755), myData5(17:1943),myData6(1:1859),myData7(1:2048),myData8(1:2048));
                        %all_data= cat(1,myData1(1:1921),myData2(1:1867),myData3(1:1682),myData4(1:1694),myData5(1:1930),myData6(1:1859),myData7(1:2048),myData8(1:2048));
                        all_data= cat(1,myData1(1:1921),myData2(1:1867),myData3(1:1682),myData4(1:1694),myData5(1:1930),myData6(1:1859),myData7(1:2048),myData8(1:2048));
                        if (size(all_data,1) ~= size(avantes_lambdas,1))
                            error('Wavelength calibration has changed, update concatenation');
                        end
                        energia = (sum(all_data) - 13E6) / 1E5;
                        energia_abs = sum(all_data); % absolute sum of values to get a reference , i.e., LA-ICP experiments in limpet 565B gave 42E6 mean value
                    end % use_avantes
                    
                    if (USE_REMOTE_PIMAX && medido_remote==1)
                        if ((now()-cuando_remote)*60*60*24 > 1) % it takes time, 2 seconds delay
                            try
                                ltemp = msrecv(sock,2); % lambdas
                                stemp = msrecv(sock,2); % spectra
                                ack = msrecv(sock,2);  % ok
                                msclose(sock);
                                medido_remote = 0;
                            catch
                                disp('Error TCP/IP accediendo a REMOTE_PIMAX - espectro');
                            end
                        end
                    end
                    
                    % refresh plot using a "fast" method  https://stackoverflow.com/questions/13102654/how-should-i-update-the-data-of-a-plot-in-matlab
                    if ~ishandle(h_disparar)
                        break;
                    end
                    if (energia > energia_limit)
                        energia_limit = energia_limit * 10;
                        subplot(9,1,7,'Parent',h_disparar);
                        xlim([0 energia_limit]);
                    end
                    subplot(9,1,7,'Parent',h_disparar);
                    h_plot_disparar.YData = all_data;
                    h_plot_disparar_barra.YData = energia;
                    yticklabels({num2str(energia_abs,'%2.1e')}); % absolute energy in exponential format
                    %set(disparar_energia_abs,'String', [ 'Abs Energy=' num2str(energia_abs,'%2.1e')]);
                    
                    %new v54 print line heights
                    subplot(9,1,8:9,'Parent',h_disparar);
                    hay_red = false;
                    hay_blue = false;
                    lineMonString='';
                    for i=1:numel(lineMonitoring_lambda)
                        %new v58 search the peaks
                        idx= lambda2pixel(avantes_lambdas,lineMonitoring_lambda(i));
                        [maxval maxidx] = max(all_data(idx-5:idx+5));
                        newidx = idx-5+maxidx-1;
                        %lineMonitoring_lambda(i) = avantes_lambdas(newidx); no need to update lambda, could drift with noise
                        val = all_data(newidx); % peak height
                        if val > ValSaturation
                            lineMonString = strcat(lineMonString,'  SAT!',lineMonitoring_name(i),':',num2str(val,'%05i ')); % magenta if value is higher than max
                            hay_red=true;
                        elseif (val>lineMonitoring_max(i))
                            lineMonString = strcat(lineMonString,'  ^',lineMonitoring_name(i),':',num2str(val,'%05i ')); % magenta if value is higher than max
                        elseif (val < 0.5*lineMonitoring_min(i))  % new v58 blue means less than 50% of min value
                            lineMonString = strcat(lineMonString,'  .',lineMonitoring_name(i),':',num2str(val,'%05i ')); % red if it is  lower than min
                            hay_blue=true;
                        else
                            lineMonString = strcat(lineMonString,'  =',lineMonitoring_name(i),':',num2str(val,'%05i ')); % just fine
                        end
                        %lineMonHeights(i) = val/((lineMonitoring_max(i)+lineMonitoring_min(i))/2)*100; % 100% = (max+min)/2
                        lineMonHeights(i) = val; % actual raw value, to better detect saturation
                    end
                    % para cambiar el color individualmente hace falta una versión más
                    % moderna de matlab:  https://es.mathworks.com/matlabcentral/answers/447017-changing-color-of-individual-bars-in-a-grouped-bar-chart
                    set(disparar_lines,'String', lineMonString); % print line heights, hopefully with colors
                    h_plot_disparar_lines.YData = lineMonHeights;
                    if hay_red
                        set(h_plot_disparar_lines,'FaceColor','r'); % v58 red means saturation
                    elseif hay_blue
                        set(h_plot_disparar_lines,'FaceColor','m'); % v58 blue changed to magenta that means less than 50% of min value
                    else
                        set(h_plot_disparar_lines,'FaceColor','g');
                    end
                    refreshdata(h_plot_disparar,'caller');
                    
                    if (USE_REMOTE_PIMAX && exist('ltemp','var')==1 && size(ltemp,1)>0)
                        %if (USE_REMOTE_PIMAX && size(ltemp,2)~=0)
                        subplot(9,1,4:6,'Parent',h_disparar);
                        %h_plot_rempimax.YData = stemp;
                        %h_plot_rempimax.XData = ltemp;
                        %refreshdata(h_plot_rempimax,'caller');
                        plot(ltemp,stemp);
                        axis([ ltemp(1) ltemp(end) 0 66000]);
                    end
                    
                    if (disparar_save_spectrum==1)
                        saved_preview_spectrum = all_data;
                        saved_preview_spectrum_lambdas = avantes_lambdas;
                        disparar_save_spectrum=0;
                    end
                    if (finishDisparar==1)
                        break
                    end
                    if (disparar_dualpulse==1)
                        DualPulseDelay = str2double(disparar_text_dualpulse.String);
                        DualPulse_Channel2_delay = 160000 - DualPulseDelay;
                        fwrite(lotis,['DELAY2=' num2str(DualPulse_Channel2_delay) char(13)]);
                        pause(0.1); if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
                        disparar_dualpulse=0;
                    end
                    if (disparar_delaycapture==1)
                        fwrite(pulsegen,['s', disparar_text_delaycapture.String, char(13)]);
                        disparar_delaycapture=0;
                    end
                    
                end % while true
                
                % close window
                if (ishandle(h_disparar)) % if closed with the red X , there is no handle
                    close(h_disparar);
                end
                if (USE_DUAL_PULSE_LASER)
                    fwrite(lotis,['QSW1OFF' char(13)]);
                    fwrite(lotis,['SHUT1OFF' char(13)]);
                    if (DualPulse) % also channel 2
                        fwrite(lotis,['QSW2OFF' char(13)]);
                        fwrite(lotis,['SHUT2OFF' char(13)]);
                    end
                    pause(0.1);
                    if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
                end
                % deactivate q-switch
                estado_qswitch=0;
                if (USE_ARDUINO_QSWITCH)  fprintf(qswitch,'0');    end
                if (USE_PULSE_GENERATOR==1) fprintf(pulsegen,['q0' char(13)]); end % activates qswitch pulses with no 'noise' pulse
                if (USE_SHUTTER_CTRL)
                    fprintf(shutter,'ens');pause(0.2);if(shutter.BytesAvailable>0);response=fread(shutter,shutter.BytesAvailable);end  % close the shutter
                end
                disparar=0;
                
            end % disparar==1
            
            if (~isempty(k) || ~ishandle(h_preview) || finishPreview == 1 || lets_go==1)  % if it has been closed
                break;
            end
            
        end % while true
        if (~ishandle(h_preview)) % if closed with the red X , there is no handle
            finishPreview = 1;
        else
            % recover paths and changes in boxes
            EXPERIMENT_DESCRIPTION = description_box.String;
            EXPERIMENT_NAME = path_box.String;
            SUBEXPERIMENT_NAME = experiment_box.String;
            experimentPath = [EXPERIMENT_NAME,'\'];
            if (USE_RECORDED_PATH ==1 && strcmp(recorded_path_filename,'')); % a specific path is not specified
                recorded_path_filename = strcat(folderPath,experimentPath,datePath,'recorded_path.mat');
            else
                recorded_path_filename = strcat(recorded_path_filename ,'\recorded_path.mat');
            end
            if (ONLY_RECORDING_PATH)
                datePath = strcat(datestr(now,'ddmmyyyy-HHMMSS'),SUBEXPERIMENT_NAME,'_PATH\');
            else
                datePath = strcat(datestr(now,'ddmmyyyy-HHMMSS'),SUBEXPERIMENT_NAME,'\');
            end
            %new v14/v64 get some measuring parameters from preview boxes
            Npulses = str2double(npulses_box.String);
            % if step is changed, change both, if not, leave it, because param2Step could be differnt and wanted to be that way
            if  param1Step ~= str2double(step_box.String)
                param1Step = str2double(step_box.String);
                param2Step = param1Step;
            end
            close(h_preview);
        end
    end % preview
    

%preview has ended
if (finishPreview==1)
    disp('Finishing...');
    %% CIERRE FINAL
    if (USE_ARDUINO_QSWITCH)
        fprintf(qswitch,'0');  % cortamos q-switch
        fclose(qswitch);
        if (USE_DUAL_PULSE_LASER==0)
            h=warndlg('¡ATENCIÓN! apagar el shutter y el Q-switch en el mando del láser.');
            uiwait(h);
        end
    end
    if (USE_LEDLIGHT)
        fclose(ledlight);
    end
    if MOTORS
        fprintf(xyz,' 1 joystick \n'); % enable joystick
        %fclose(theta);
        fclose(xyz);
    end
    if (USE_DUAL_PULSE_LASER)
        fwrite(lotis,['QSW1OFF' char(13)]);
        fwrite(lotis,['SHUT1OFF' char(13)]);
        if (DualPulse) % also channel 2
            fwrite(lotis,['QSW2OFF' char(13)]);
            fwrite(lotis,['SHUT2OFF' char(13)]);
        end
        if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
        pause(0.1);
        fwrite(lotis,['STOP' char(13)]);    % new v55 stop the lamp, so the experiments could be left running
        pause(0.1);
        if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
        fclose(lotis);
    end
    if (USE_PULSE_GENERATOR==1) fprintf(pulsegen,['q0' char(13)]); fclose(pulsegen); end % activates qswitch pulses with no 'noise' pulse
    if (USE_SHUTTER_CTRL)
        fprintf(shutter,'ens?');pause(0.1);if(shutter.BytesAvailable>0);response=fread(shutter,shutter.BytesAvailable);end
        if (numel(response) >= 6 && response(6)=='1')
            fprintf(shutter,'ens');pause(0.1);if(shutter.BytesAvailable>0);response=fread(shutter,shutter.BytesAvailable);end % open the shutter
        end
        fclose(shutter);
    end
    
    stop(vid);
    if (TAKE_WEBCAM_PHOTO==1)
        stop(vid2);
    end
    if (havePolygon || haveLinePath)  || (exist('h_stitch','var')  && ~StitchedIsLoadedNotCreated)
        sino= questdlg(['There is a StitchedImage/LinePath/polygon created, do you want to save it?' ], 'Yes', 'No');
        switch sino
            case 'Yes'
                % creates folder for experiment
                mkdir(strcat(folderPath,experimentPath,datePath));
                if havePolygon
                    disp(['Saving polygon to ' strcat(folderPath,experimentPath,datePath,'polygon.mat')]);
                    save (strcat(folderPath,experimentPath,datePath,'polygon.mat'),'roiPol');
                end
                if haveLinePath
                    save (strcat(folderPath,experimentPath,datePath,'linePath.mat'),'roiLine');
                    disp(['Saving linePath to ' strcat(folderPath,experimentPath,datePath,'linePath.mat')]);
                end
                if ishandle(h_stitch)
                    savefig(h_stitch,strcat(folderPath,experimentPath,datePath,'stitchedImage.fig'));
                    imwrite(stitchImage, strcat(folderPath,experimentPath,datePath,'stitchedImage.png')); % full resolution image
                end
            case 'No'
                % do nothing
        end
    end
    try
        close(h_danger); % danger windows uses to hang around
    catch
    end
    disp('Done!');
    
    return; % actually, it ends
end

%%%%%%%% aquí se empieza a medir!!

%% get current coordinates 
    try
        fprintf(xyz,' p \r'); pause(0.05);
        current_coords = fscanf(xyz,'%f %f %f'); % this reads the current coordinates x y z
        backToStartingCoordinates = current_coords; % new v62 at the end, goes back here
    catch
    end

timePerPoint = (Npulses *0.1) + 2.33 ; % the npulsesx0.1secs plus "overhead", redefined here just in case
%newv61 if havePolygon, change 2D corners to fit the shape
if havePolygon && absoluteMovement
    xmin=min(roiPol.Position(:,1));
    xmax=max(roiPol.Position(:,1));
    ymin=min(roiPol.Position(:,2));
    ymax=max(roiPol.Position(:,2));
    if strcmp(param1Name,'posX') param1From = xmin; param1To = xmax;  end
    if strcmp(param2Name,'posX') param2From = xmin; param2To = xmax;  end
    if strcmp(param1Name,'posY') param1From = ymin; param1To = ymax;  end
    if strcmp(param2Name,'posY') param2From = ymin; param2To = ymax;  end 
    %new v14/v64 calcs again effecti number of points
    xmin=min(roiPol.Position(:,1));
    xmax=max(roiPol.Position(:,1));
    ymin=min(roiPol.Position(:,2));
    ymax=max(roiPol.Position(:,2));
    ejex = xmin:param1Step:xmax;
    ejey = ymin:param2Step:ymax;
    [gridx,gridy]=meshgrid(ejex,ejey);
    effectivePoints =sum(inROI(roiPol,gridx,gridy));
    disp('Warning!! 2D measurement corners changed to polygon shape');
    disp(['  Estimated time: ' datestr(seconds(timePerPoint*effectivePoints),'DD:HH:MM:SS')]);
end
if haveLinePath && ~absoluteMovement % a measurement path has been defined but last loaded/created was not a polygon
     disp(['Calculating the measurement path using deltaR = ' num2str(param1Step) ' ... ']);
     deltaR = param1Step;
     USE_RECORDED_PATH=1;
     recorded_path = zeros(0,5); % rebuild the relative movements
     rp_realPos = zeros(0,2); % the same vector as recorded_path but in real position, so we can draw point on the stitched image
     rp_AbsStartingPoint= zeros(2,1);
     rp_AbsStartingPoint(1) = roiLine.Position(1,1); % new v09 this is absolute starting point to move to under use_recorded_path && haveLinePath
     rp_AbsStartingPoint(2) = roiLine.Position(1,2);
     fprintf(xyz,' p \r'); pause(0.05);
     current_coords = fscanf(xyz,'%f %f %f'); % this reads the current coordinates x y z
     if ishandle(h_stitch)
         figure(h_stitch);
     end
     nTotalPointsDeltaR=0;
     for i=1:size(roiLine.Position,1)-1
         nextX = roiLine.Position(i+1,1);
         nextY = roiLine.Position(i+1,2);
         if size(rp_realPos,1) > 0
            prevX = rp_realPos(end,1); % the last point could no be the same due to deltaR
            prevY = rp_realPos(end,2);
         else
            prevX = rp_AbsStartingPoint(1); % for the first point
            prevY = rp_AbsStartingPoint(2);
         end
         if nextX>prevX  dirX=1; else dirX=-1; end
         if nextY>prevY  dirY=-1;else dirY=1;  end
         mantrack_dist =  sqrt((nextX-prevX)^2+(nextY-prevY)^2); % distance in mm
         mantrack_arctan = atan(abs(nextY-prevY)/abs(nextX-prevX));
         npoints = floor(mantrack_dist/ deltaR);
         for j=1:npoints
            recorded_path = [ recorded_path ; -deltaR*sin(mantrack_arctan)*dirY, deltaR*cos(mantrack_arctan)*dirX, 0, 0 , 0];
            prevX = prevX + deltaR*cos(mantrack_arctan)*dirX;
            prevY = prevY - deltaR*sin(mantrack_arctan)*dirY;
            rp_realPos = [ rp_realPos ; prevX  , prevY ];
            if ishandle(h_stitch)
                drawpoint('Position',[rp_realPos(end,1),rp_realPos(end,2)],'Color','yellow','MarkerSize',3);
            end
            nTotalPointsDeltaR = nTotalPointsDeltaR +1;
         end
         drawnow;
     end % for i
     old_recorded_path = recorded_path;
     if size(recorded_path,1) > size(param1From:param1Step:param1To,2) % the number of points is larger than the allocated buffer
         disp('Warning!! The length of the deltaR loop is smaller than the number of recorded points, param1To has been changed');
         param1To = param1From+ deltaR*(size(recorded_path,1)+1);
     end
     effectivePoints = nTotalPointsDeltaR;
end % if haveLinePath

Nparam1 = size(param1From:param1Step:param1To,2); % OJO! v55 changed from round((param1To-param1From)/param1Step +1 );
Nparam2 = size(param2From:param2Step:param2To,2); % round((param2To-param2From)/param2Step +1 ) ;
Npoints = Nparam1 * Nparam2 ;
point = 1 ; % index for incremental file naming
pointCounterEffective = 1; % index for effective points (measured, not skipped)
param1Index = 1; param2Index = 1;
ratio_MgCa_toplot = NaN(Nparam1,Nparam2); % se calcula el ratio espectro a espectro
SNR_Mg = NaN(Nparam1,Nparam2); % para el cálculo de la SNR
SNR_Ca = NaN(Nparam1,Nparam2);
satlines_pix = zeros(size(satlines,2),1);
checkSurface_matrix = zeros(Nparam1, Nparam2); %new v55 store 1 for spatial points actually measured
spectrometerTemperature = zeros(Nparam1, Nparam2); % new v57 to store spectrometer temperature (internal NTC in Avantes)
effectivePoints = Nparam1*Nparam2; % default number of points for absolute measurement       

whenStarted = now;

%new v62 if NO_SHOW_IMAGE=0 && Npoints > 2000 -> problema de memory leak
if NO_SHOW_IMAGE==0 && Npoints > 2000 && 0 % % disabled!! becasue the memory leak bug is fixed
    sino= questdlg(['Warning!! number of points is > 2000, image in monitor window WILL BE DISABLED, photo files (if configured) are stored anyway, proceed?' ], 'Yes', 'No');
    switch sino
        case 'Yes'
            NO_SHOW_IMAGE=1;
        case 'No'
            % do nothing
    end
end


% creates folder for experiment
mkdir(strcat(folderPath,experimentPath,datePath));
displaysave([EXPERIMENT_DESCRIPTION,'\n']);
try % new v65 mfilename is empty and we don't know why
    copyfile([ mfilename('fullpath') '.m'] ,[ folderPath experimentPath datePath mfilename '.m']); %copy the current program to the experiment folder, easier to repeat
catch
    disp('Warning! mfilename() gave error, the script has not been copied to the experiment folder.');
end

% if USE_RECORDED_PATH the variable is ...

if (USE_LOCAL_PIMAX)
    % PULSES per experiment
    objExp.SetParam('EXP_FILEINCENABLE',0); % INCREMENT COUNT SHOULD BE DISABLE FOR THE SAVE COMMAND TO WORK
    if (AUTOCOLLECT==1)
        objExp.SetParam('EXP_SEQUENTS',1); % only one spectrum per experiment
    else
        objExp.SetParam('EXP_SEQUENTS',Npulses);
    end
    % gain
    actual_gainIntensifier = gainIntensifier;
    previous_gain = gainIntensifier; % incremental gain changes due to autogain
    objExp.SetParam('EXP_INTENSIFER_GAIN',actual_gainIntensifier);
    
    % pulser
    objPul.SetParam('TGC_PULSE_DELAY',delayCapture);
    objPul.SetParam('TGC_TRIG_FREQUENCY',frequencyPulses);
    objPul.SetParam('TGC_PULSE_WIDTH',gateWidth);
    objPul.Process('TGP_SET_GATE_PULSE'); % commits the software parameters into hardware
    
    % get background and save
    if (BACKGROUND==1)
        displaysave('Getting background ...  close the lase shutter and press any key...');
        pause;
        objExp.SetParam('EXP_DARKNAME',strcat(folderPath,experimentPath,datePath,backgroundFile));
        %OJO!!! background does not shut off the laser
        if (USE_ARDUINO_QSWITCH)  fprintf(qswitch,'0');    end
        objExp.AcquireBackground; %% this command also sets backgroundsubstraction = true
        if (USE_ARDUINO_QSWITCH)  fprintf(qswitch,'0');    end
        objExp.SetParam('EXP_BBACKSUBTRACT','1');  % parameter not documented: change data_correction setting in experiment setup
        % restore correct delay
        objPul.SetParam('TGC_PULSE_DELAY',delayCapture);
        disp('Press any key to start...');
        pause;
    else
        objExp.SetParam('EXP_BBACKSUBTRACT','0');  % parameter not documented: change data_correction setting in experiment setup
    end
    
    % save experiment setup parameter just in case
    displaysave('Saving experiment setup');
    objExp.SetParam('EXP_SETUPNAME',strcat(folderPath,experimentPath,datePath,setupFile));
    objExp.Save;
end % PIMAX

if(USE_REMOTE_PIMAX==1)
    try
        % npulses could be 1 if preview-disparar was used
        sock = msconnect(remote_pimax_ip,3000);
        comando = [ 'npulses=' num2str(remote_pimax_npulses)];
        ret=mssend(sock, comando, 2);ack = msrecv(sock,3); msclose(sock);
        if (strcmp(ack,'ok')==0)
            warndlg('¡ATENCIÓN! no hay respuesta del servidor REMOTE_PIMAX');
        end
    catch
        disp('Error TCP/IP accediendo a la cámara PIMAX3 remota.');
    end
end


scrsz=get(0,'ScreenSize');
h_focus = figure('Position',[scrsz(3)/2 40 scrsz(3)/2 scrsz(4)-100],'Name','Measure monitoring');
esc_button = uicontrol('Style', 'PushButton', 'String', 'Finish', 'Callback', 'finish=1;', 'Position',[10 10 50 20],'ForegroundColor','red');
pau_button = uicontrol('Style', 'PushButton', 'String', 'Pause', 'Callback', 'pausing=1;', 'Position',[100 10 50 20]);
skip_button = uicontrol('Style', 'PushButton', 'String', 'Skip', 'Callback', 'skipping=1;', 'Position',[200 10 50 20]);
autofocus_button = uicontrol('Style', 'PushButton', 'String', 'Autofocus', 'Callback', 'hacer_autofocus=1;', 'Position',[300 10 80 20]);
save_button = uicontrol('Style', 'PushButton', 'String', 'Save data', 'Callback', 'do_saveData=1;', 'Position',[400 10 80 20]);

hacer_autofocus=0;
defocus_mm = 0;
do_saveData=0;


% reference frame for autofocus
if (CAMERA)
    % new v28 imshow ONLY ONCE at the beginning
    % changed v38 the problem was a memory leak already solved
    frame = captureFrameGuppy; % get new image
    %figure(h_focus);
    set(0,'CurrentFigure',h_focus);
    subplot(5,1,1);
    h_ims_afocus=imshow(frame,'InitialMagnification', 67);
    
    if (AUTOFOCUS && strcmp(AUTOFOCUS_METHOD,'image')==1) % nothing to do because this method does not uses reference focus
    end
    
    if (AUTOFOCUS && strcmp(AUTOFOCUS_METHOD,'laserspot')==1)  % get the reference
        ref_pixel = GetXYfocus(theta,vid,pixel_shotX,pixel_shotY);
        if (ref_pixel == [0 0])
            error('reference spot for focus not found, move sample in Z to see the laser spot');
        else
            displaysave(['Reference focus X=' int2str(ref_pixel(1)) ' ' 'Y=' int2str(ref_pixel(2)) ]);
        end
    end
    
    % guardamos la imagen antes de empezar
    if (~ ONLY_RECORDING_PATH) % if only_recording_path no photos to speed up things
        frameBefore = captureFrameGuppy;
        imwrite(frameBefore,strcat(folderPath,experimentPath,datePath,'BeforePhoto.png')); %save current frame with the same name as the experiment
    end
end

% webcam image
if (TAKE_WEBCAM_PHOTO)
    if (~ ONLY_RECORDING_PATH) % if only_recording_path, no photos to speed up things
        frame = captureFrameWebcam;
        imwrite(frame,strcat(folderPath,experimentPath,datePath,filename,'_webcam.png')); %save current frame with the same name as the experiment
    end
end

%%the code for SCANLINE==2 goes here
if (SCANLINE==2)
    sound(sound2);sound(sound3); % sound alert
    h_rubber=figure('Name','Current image - draw first SCANLINE over the surface'); % create new large window
    ok_rubber_button = uicontrol('Style', 'PushButton', 'String', 'Continuar', 'Callback', 'salir_rub=1;', 'Position',[10 10 50 20]);
    %again_rubber_button = uicontrol('Style', 'PushButton', 'String', 'Otra vez', 'Callback', 'again_rub=1;', 'Position',[200 10 50 20]);
    length_rubber_text = uicontrol('Style', 'text', 'String', 'Length (mm) =', 'Position',[100 10 200 20]);
    angle_rubber_text = uicontrol('Style', 'Text', 'String', 'Angle º(º) =', 'Position',[300 10 200 20]);
    
    frame = captureFrameGuppy;
    h_ims_rub=imshow(frame,'InitialMagnification', 67); hold on;
    target(pixel_shotX,pixel_shotY); % draw a cross around
    if (ESCALA==1) % draw 1mm ticks starting at the laser spot
        if calibrated_pixels_mm_XY_X > 300   % new v58
            mm1o01=calibrated_pixels_mm_XY_X*0.1; % line every 100um
            numberOflines=7;
        else
            mm1o01=calibrated_pixels_mm_XY_X; % line every 1mm
            numberOflines=3;
        end
        hold on;
        for horticks=1:numberOflines % cortos hacia la izquierda
            line([pixel_shotX-(horticks-1)*mm1o01  pixel_shotX-(horticks-1)*mm1o01],[0 300],'Color','w');
            line([pixel_shotX-(horticks-1)*mm1o01+1  pixel_shotX-(horticks-1)*mm1o01+1],[0 300],'Color','k');
        end
        for horticks=1:numberOflines  % largos hacia la derecha
            line([pixel_shotX+(horticks)*mm1o01  pixel_shotX+(horticks)*mm1o01],[0 600],'Color','w');
            line([pixel_shotX+(horticks)*mm1o01+1  pixel_shotX+(horticks)*mm1o01+1],[0 600],'Color','k');
        end
    end
    salir_rub=0;
    h_rub_l = line ([0 0] , [0 0]);
    h_rub = impoint();
    while(salir_rub==0)
        frame = captureFrameGuppy;
        set(h_ims_rub,'CData',frame); % new! v28 update image without a new imshow, promise speed things up
        target(pixel_shotX,pixel_shotY); % draw a cross around
        if (ESCALA==1) % draw 1mm ticks starting at the laser spot
            if calibrated_pixels_mm_XY_X > 300   % new v58
                mm1o01=calibrated_pixels_mm_XY_X*0.1; % line every 100um
                numberOflines=7;
            else
                mm1o01=calibrated_pixels_mm_XY_X; % line every 1mm
                numberOflines=3;
            end
            hold on;
            for horticks=1:numberOflines % cortos hacia la izquierda
                line([pixel_shotX-(horticks-1)*mm1o01  pixel_shotX-(horticks-1)*mm1o01],[0 300],'Color','w');
                line([pixel_shotX-(horticks-1)*mm1o01+1  pixel_shotX-(horticks-1)*mm1o01+1],[0 300],'Color','k');
            end
            for horticks=1:numberOflines  % largos hacia la derecha
                line([pixel_shotX+(horticks)*mm1o01  pixel_shotX+(horticks)*mm1o01],[0 600],'Color','w');
                line([pixel_shotX+(horticks)*mm1o01+1  pixel_shotX+(horticks)*mm1o01+1],[0 600],'Color','k');
            end
        end
        pos = getPosition(h_rub);
        rub_lengthx = (pos(1)-pixel_shotX);
        rub_lengthy = (pos(2)-pixel_shotY);
        rub_length = 2* sqrt((pos(1)-pixel_shotX)^2+(pos(2)-pixel_shotY)^2); % [x y]
        rub_length_mm = 2*sqrt((rub_lengthx/calibrated_pixels_mm_XY_X)^2+(rub_lengthy/calibrated_pixels_mm_XY_Y)^2);
        rub_angle = atand ((pixel_shotY-pos(2))/(pos(1)-pixel_shotX));
        %rub_length = sqrt((pos(1,1)-pos(2,1))^2+(pos(1,2)-pos(2,2))^2); % pos(1,1)=x1 pos(2,1)=x2
        %rub_angle = atand ((pos(1,2)-pos(2,2))/(pos(2,1)-pos(1,1)));
        set(length_rubber_text,'String', [ 'Length(mm)=' num2str(rub_length_mm,'%3.2f')]);
        set(angle_rubber_text,'String', [  'Angle(º)=' num2str(rub_angle,'%3.1f')]);
        x1 = pos(1);
        x2 = pixel_shotX - (pos(1)-pixel_shotX);
        y1 = pos(2);
        y2 = pixel_shotY + (pixel_shotY-pos(2));
        % setPosition(h_rub,[x1 y1; x2 y2]);
        delete(h_rub_l);
        h_rub_l = line ([x1 x2] , [y1 y2]);
        drawnow(); % needed to process the event qeue
        if (isempty(h_rubber) || ~ishandle(h_rubber) || salir_rub==1)  % if it has been closed
            break;
        end
    end
    
    close(h_rubber);
    scanLineLength = rub_length_mm;
    scanLineAngle = rub_angle;
    displaysave(['SCANLINE=2  length(mm)=' num2str(rub_length_mm,'%3.2f') ' angle(º)=' num2str(rub_angle,'%3.1f') ]);
    
end

% this is the "master" array storing the raw captured spectra. if > 300MB we'll try to split param2 in individual files
if (USE_AVANTES)
    Npixeles = Npixeles_avantes;
else
    Npixeles = Npixeles_pimax;
end
if ((Npixeles*Npulses*Nparam1*Nparam2 > 600E6/4 && Nparam2 > 1) && SPLIT_PARAM2==1) % v55 changed to AND, so real limit can be tested
    h=warndlg('¡ATENCIÓN! El número de datos es muy elevado, se van a guardar con SPLIT_PARAM2=1');
    SPLIT_PARAM2 = 1;
    if (MULTI_WINDOW)
        MW_lambdas = zeros(Npixeles,size(MW_centralwl,2)); % realiza varias medidas a las lambdas indicadas en MW_centralwl. Guarda los espectros en MW_lambdas[:,num_ventanas] y MW_spectra(:,npulses,idx1,idx2,nw)
        MW_spectra =  zeros(Npixeles,Npulses, Nparam1, 1,size(MW_centralwl,2),'single');
    else
        spectra = zeros(Npixeles,Npulses, Nparam1, 1,'single');
    end
    if (USE_OCEANOPTICS==1)
        hr2000_spectra = zeros(Npixeles_hr2000,Npulses, Nparam1, 1,'single');
    end
    %avantes_spectra = zeros(Npixeles_avantes,Npulses, Nparam1,
    %Nparam2,'single'); variable removed as it takes a lot of memory and it
    %is redundant
else
    if (MULTI_WINDOW)
        MW_lambdas = zeros(Npixeles,size(MW_centralwl,2)); % realiza varias medidas a las lambdas indicadas en MW_centralwl. Guarda los espectros en MW_lambdas[:,num_ventanas] y MW_spectra(:,npulses,idx1,idx2,nw)
        MW_spectra =  zeros(Npixeles,Npulses, Nparam1, Nparam2,size(MW_centralwl,2),'single');
    else
        %new v65 gave errors for size > 31GB, because was limited to 100% of physical RAM (default setting preferences->workspcae)
        try
          spectra = zeros(Npixeles,Npulses, Nparam1, Nparam2,'single'); % spectra() es el que se usa para RSD, ratios sobre la marcha ...
        catch
          error('spectra() size > physical RAM and gave and cannot continue, change preferences->workspace');  
        end
        %new v58 spectra8raw stores the complete CCD spectra
        if STORE_SPECTRA8RAW
            spectra8raw = zeros(2048,Npulses,Nparam1,Nparam2,8,'single');
        else
            if exist('spectra8raw','var')
                clear('spectra8raw'); % delete it, if the problem was too much memory before, we need to clear it
            end
        end
        if (USE_OCEANOPTICS==1)
            hr2000_spectra = zeros(Npixeles_hr2000,Npulses, Nparam1, Nparam2,'single');
        end
        %avantes_spectra = zeros(Npixeles_avantes,Npulses, Nparam1, Nparam2,'single');
    end
end
if (USE_REMOTE_PIMAX==1)
    spectra_remote_pimax = zeros(Npixeles_pimax, remote_pimax_npulses, Nparam1, Nparam2,'single');
    lambdas_remote_pimax = zeros(Npixeles_pimax,1);
end

spectra_is_valid = zeros(Npulses, Nparam1, Nparam2);  % =1 for the spectra labeled as valids
num_valid_spectra_matrix = zeros(Nparam1, Nparam2);  % each value is the number of valid spectra for that experiment
autoGains = zeros(Nparam1, Nparam2);  % chaging autogains
lambdas = zeros(Npixeles,1);
if USE_AVANTES
    lambdas = avantes_lambdas;
end

% ... and processing results
% global variables for processing
Rsd_1pixel = zeros(Nparam1, Nparam2,1);
Rsd_ratio  = zeros(Nparam1, Nparam2,1);
Rsd_1pixel_mean = zeros(Nparam1, Nparam2,1);
Rsd_ratio_mean = zeros(Nparam1, Nparam2,1);
all_intensities = zeros(Npulses,Nparam1,Nparam2,1);
pixel_captured_from_first_spectrum = 0; % the first valid spectrum will be used to calculate the pixels for 1pixel and ratio RSD

axisParam1 = zeros(Nparam1,1);  axisParam2 = zeros(Nparam2,1); % plot axis
first_time_lambdas = 1; % read the lambda vector only once
spatial_points = zeros(Nparam1,Nparam2,3); % store the spatial X,Y,Z of every point, so the results can be spatially mapped
spatial_points(1,1,1)=0;spatial_points(1,1,2)=0;spatial_points(1,1,3)=0;  % first point is (0,0,0) in spatial coordinates
timer_experiment=tic; % the first time we need a time reference
now1=toc(timer_experiment); now2=now1; now3=now1; % we need a value for this vars 1st time
skipThisPoint = 0; % with ROI, some points are not measured in 2D configuration
panicAvantes=0; % set to 1 if ready takes too much time, it hangs up sometimes new v62
counterAvantesErrors=0; % new v64 retry several times
%new v59 fscanf get a timeout error once in a while, could be a memory problem or a driver issue, flowcontrol is set to 'software' and I don't remember why
%new v61 position is read once at the beginning
if MOTORS
    try
        fprintf(xyz,' p \r');pause(0.05);
        current_coords = fscanf(xyz,'%f %f %f'); % this reads the current coordinates x y z
        absX = current_coords(2);
        absY = current_coords(1);
        absZ = current_coords(3);
    catch
    end
    %new v62
    % in MANUAL_TRACKING && USE_RECORDED_PATH mode, coordinates are relative, we should move first to the starting point of the pre-recordedpath
    if USE_RECORDED_PATH && haveLinePath % in this case, we have an absolute starting point to move to
        moveString = sprintf(' %0.3f %0.3f %0.3f m \r',rp_AbsStartingPoint(2), rp_AbsStartingPoint(1), absZ); % absolute values, command move (m)
        fprintf(xyz,moveString); 
        justAbsMovDone=tic(); % to check if motors are effectively stopped
        statusXYZ=1; % v64 code to wait for motors to stop, already tested
        while statusXYZ==1        %is moving,
            pause(0.05);
            if(xyz.BytesAvailable>0);fread(xyz,xyz.BytesAvailable);end
            fprintf(xyz,' status \r');
            statusXYZ = fscanf(xyz,'%f'); % status word
        end
    end
end % if MOTORS


%%%%%%%%%%%%%%%%%% MASTER LOOPS 
for param2IndexReal= 1:Nparam2
    
    if (SPLIT_PARAM2==1)
        param2Index = 1;
        spectra(:,:,:,1) = 0; % split mode: before next param2, clear spectra
    else
        param2Index = param2IndexReal;
    end
    
    checkSurfaceAlreadyPositive = false; % let's start with no measurment
    for param1Index= 1:Nparam1
        param1 = param1From + (param1Index-1)*param1Step; axisParam1(param1Index) = param1;
        param2 = param2From + (param2IndexReal-1)*param2Step; axisParam2(param2IndexReal) = param2;
        
        %% CHANGE THINGS ACCORDING TO THE NAME OF PARAM1 & PARAM2 new! v47
        
        % 'deltaT' for every spatial point, perform a transverse movement
        if ( strcmp(param1Name,'deltaT') == 1 && param1Index == 1)
            dx2D =  L2D/2 * cos(pi/180*DeltaT_ANGLE);
            dy2D =  L2D/2 * sin(pi/180*DeltaT_ANGLE);
            displaysave(['deltaT first movement XY ' num2str(dx2D) ' ' num2str(dy2D)]);
            moveString = sprintf(' %0.3f %0.3f %0.3f r \r',dy2D, dx2D, 0); %
            fprintf(xyz,moveString); %move axis relative
        end
        % 'DualPulseDelay_ns' change interpulse delay of lotis laser
        if ( strcmp(param1Name,'DualPulseDelay_ns') == 1)
            DualPulseDelay = param1;
            DualPulse_Channel2_delay = 160000 - DualPulseDelay; % fixed for optimum lamp charge
            fwrite(lotis,['DELAY2=' num2str(DualPulse_Channel2_delay) char(13)]);
            pause(0.1); if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
            displaysave(['LOOP CHANGE idx1: DualPulseDelay_ns=' num2str(DualPulseDelay)]);
        end
        % 'DualPulseDelay_ns' change interpulse delay of lotis laser
        if ( strcmp(param1Name,'DualPulseDelay_ns') == 2)
            DualPulseDelay = param2;
            DualPulse_Channel2_delay = 160000 - DualPulseDelay; % fixed for optimum lamp charge
            fwrite(lotis,['DELAY2=' num2str(DualPulse_Channel2_delay) char(13)]);
            pause(0.1); if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
            displaysave(['LOOP CHANGE idx2: DualPulseDelay_ns=' num2str(DualPulseDelay)]);
        end
        % 'delayCaptureForReal_ns' changes the effective capture delay
        if ( strcmp(param1Name,'delayCaptureForReal_ns') == 1)
            delayCaptureForReal = param1;
            if (USE_DUAL_PULSE_LASER)
                delayCaptureAvantes = 160400 - 1300 + delayCaptureForReal; % [ns]
            else
                delayCaptureAvantes = 160400 - 1300 + delayCaptureForReal; % [ns]
            end
            % write the changes to pulsegen for avantes
            fwrite(pulsegen,['s' num2str(delayCaptureAvantes) char(13)]);   % delay of spectrometer
            pause(0.1); if(pulsegen.BytesAvailable>0);fread(pulsegen,pulsegen.BytesAvailable);end
            remote_pimax_delayCapture = delayCaptureAvantes + 1300; % el avantes tarda realmente 1300ns en capturar desde el trigger
            %write to pimax , should be checked that it works because i am not sure some commands are effe
            if (USE_REMOTE_PIMAX==1)
                try
                    sock = msconnect(remote_pimax_ip,3000);
                    comando = [ 'delay=' num2str(remote_pimax_delayCapture)];
                    ret=mssend(sock, comando, 2);ack = msrecv(sock,3); msclose(sock);
                    if (strcmp(ack,'ok')==0)
                        warndlg('¡ATENCIÓN! no hay respuesta del servidor REMOTE_PIMAX');
                    end
                catch
                    disp('Error TCP/IP accediendo a la cámara PIMAX3 remota.');
                end
            end
            displaysave(['LOOP CHANGE idx1: delayCaptureForReal_ns=' num2str(delayCaptureForReal)]);
        end
        % 'delayCaptureForReal_ns' changes the effective capture delay
        if ( strcmp(param2Name,'delayCaptureForReal_ns') == 1)
            delayCaptureForReal = param2;
            if (USE_DUAL_PULSE_LASER)
                delayCaptureAvantes = 160400 - 1300 + delayCaptureForReal; % [ns]
            else
                delayCaptureAvantes = 160400 - 1300 + delayCaptureForReal; % [ns]
            end
            % write the changes to pulsegen for avantes
            fwrite(pulsegen,['s' num2str(delayCaptureAvantes) char(13)]);   % delay of spectrometer
            pause(0.1); if(pulsegen.BytesAvailable>0);fread(pulsegen,pulsegen.BytesAvailable);end
            remote_pimax_delayCapture = delayCaptureAvantes + 1300; % el avantes tarda realmente 1300ns en capturar desde el trigger
            %write to pimax , should be checked that it works because i am not sure some commands are effe
            if (USE_REMOTE_PIMAX==1)
                try
                    sock = msconnect(remote_pimax_ip,3000);
                    comando = [ 'delay=' num2str(remote_pimax_delayCapture)];
                    ret=mssend(sock, comando, 2);ack = msrecv(sock,3); msclose(sock);
                    if (strcmp(ack,'ok')==0)
                        warndlg('¡ATENCIÓN! no hay respuesta del servidor REMOTE_PIMAX');
                    end
                catch
                    disp('Error TCP/IP accediendo a la cámara PIMAX3 remota.');
                end
            end
            displaysave(['LOOP CHANGE idx2: delayCaptureForReal_ns=' num2str(delayCaptureForReal)]);
        end
        % 'gain' changes the gain of PIMAX camera
        if ( strcmp(param1Name,'gain') == 1)
            gainIntensifier = param1;
            remote_pimax_gainIntensifier = gainIntensifier;
            %write to pimax , should be checked that it works because i am not sure some commands are effe
            if (USE_REMOTE_PIMAX==1)
                try
                    sock = msconnect(remote_pimax_ip,3000);
                    comando = [ 'gain=' num2str(remote_pimax_gainIntensifier)];
                    ret=mssend(sock, comando, 2);ack = msrecv(sock,3); msclose(sock);
                    if (strcmp(ack,'ok')==0)
                        warndlg('¡ATENCIÓN! no hay respuesta del servidor REMOTE_PIMAX');
                    end
                catch
                    disp('Error TCP/IP accediendo a la cámara PIMAX3 remota.');
                end
            end
            displaysave(['LOOP CHANGE idx1: gain=' num2str(gainIntensifier)]);
        end
        % 'gain' changes the gain of PIMAX camera
        if ( strcmp(param2Name,'gain') == 1)
            gainIntensifier = param2;
            remote_pimax_gainIntensifier = gainIntensifier;
            %write to pimax , should be checked that it works because i am not sure some commands are effe
            if (USE_REMOTE_PIMAX==1)
                try
                    sock = msconnect(remote_pimax_ip,3000);
                    comando = [ 'gain=' num2str(remote_pimax_gainIntensifier)];
                    ret=mssend(sock, comando, 2);ack = msrecv(sock,3); msclose(sock);
                    if (strcmp(ack,'ok')==0)
                        warndlg('¡ATENCIÓN! no hay respuesta del servidor REMOTE_PIMAX');
                    end
                catch
                    disp('Error TCP/IP accediendo a la cámara PIMAX3 remota.');
                end
            end
            displaysave(['LOOP CHANGE idx2: gain=' num2str(gainIntensifier)]);
        end
        if absoluteMovement == 1 % using posX, posY, posZ, movement should be done BEFORE measurement
            if strcmp(param1Name,'posX')
                absX = param1; % go-to position is current idx1 value
            end
            if strcmp(param2Name,'posX')
                absX = param2;
            end
            if strcmp(param1Name,'posY')
                absY = param1;
            end
            if strcmp(param2Name,'posY')
                absY = param2;
            end
            if strcmp(param1Name,'posZ')
                absZ = param1;
            end
            if strcmp(param2Name,'posZ')
                absZ = param2;
            end
            % we should move here, BUT if havePolygon==1, we need to check ROI before
            % new v65 skiping can be due to inROI or skipFunction below
 
            if isequal(skipFunction,'roi') && havePolygon &&  ~inROI(roiPol,absX,absY)
                displaysave(['Skipping point outside ROI ' num2str(absX) ',' num2str(absY)]);
                skipThisPoint=1;
            elseif (~isequal(skipFunction,'') && ~isequal(skipFunction,'none')) && ~isequal(skipFunction,'roi')
                skipThisPoint = skipFunctionHandler(absX,absY,param1From,param1To,param1Step,param2From,param2To,param2Step);
                if skipThisPoint
                    displaysave(['Skipping point due to skipFunction ' skipFunction ' at ' num2str(absX) ',' num2str(absY)])
                end
            else
                skipThisPoint = 0; % do the thing
            end
            if ~skipThisPoint
                moveString = sprintf(' %0.3f %0.3f %0.3f m \r',absY, absX, absZ); % absolute values, command move (m)
                fprintf(xyz,moveString); %move axis
                justAbsMovDone=tic(); % to check if motors are effectively stopped
                skipThisPoint=0;
                statusXYZ=1; % v64 code to wait for motors to stop, already tested
                while statusXYZ==1        %is moving,
                    pause(0.05);
                    if(xyz.BytesAvailable>0);fread(xyz,xyz.BytesAvailable);end
                    fprintf(xyz,' status \r');
                    statusXYZ = fscanf(xyz,'%f'); % status word
                end
                
            end
        end
        
        
        % set datafile , mar2017 add multiwindow option
        filename = strcat(sprintf('%04.0f',point),'_',param1Name,'=',sprintf('%0.2f',param1),'_',param2Name,'=',sprintf('%0.2f',param2));
        filename = strrep(filename,'.',','); %replaces dots with commas, just in case
        displaysave(filename);
        
        if (USE_LOCAL_PIMAX)
            objExp.SetParam('EXP_DATFILENAME',strcat(folderPath,experimentPath,datePath,filename,'.spe'));
            if (AUTOCOLLECT)  % only in the autocollect mode the gain is reset to the original value for every spot
                actual_gainIntensifier = gainIntensifier;
                objExp.SetParam('EXP_INTENSIFER_GAIN',gainIntensifier); % every new point the gain is restored
            end
        end
        if(finish)
            displaysave('Breaking from after-experiment');
            break
        end
        
        % in the 2D mode (param1Name = 'deltaT'), we need to move to one
        % extreme BEFORE the first point
        
        
        if(pausing)
            fprintf(xyz,' 1 joystick \n');
            displaysave('Pause ... press button again to continue, laser armed?');
            pausing=0;
            set(pau_button , 'string','un-pause');
            while (pausing==0)
                drawnow;
            end
            set(pau_button , 'string','Pause');
            pausing = 0;
        end
        if do_saveData % new v13/64 save data so far, to peek what is going on
            displaysave('Saving data due to user request... it can take a while...');
            save (strcat(folderPath,experimentPath,datePath,dataFile),'-regexp','^(?!(wrapper|vid*|obj*|stitchImage)$).','-v7.3');
            do_saveData=0;
        end
        if ~skipThisPoint
            pause(SLEEP_BETWEEN_CAPTURES) % new v61 I believe that this pause should be here, after absolute movement.
        end
        % perform autofocus if requested
        if (CAMERA && ~skipThisPoint)
            %figure(h_focus);
            set(0,'CurrentFigure',h_focus);
            
            if (AUTOFOCUS && strcmp(AUTOFOCUS_METHOD,'image')==1) % performs an initial autofocus
                if (SOUND_STOP==0 || (SOUND_STOP==1 && param1Index==1)) % in SOUND_STOP mode, the focus should be done only once in every param2 point
                    focusLastZ = 0; % start value for z
                    focus_antes = minimizaFocus(0);
                    options = optimset('Display','notify','TolX',0.001,'MaxFunEvals',20);
                    [z focusVal, exitFlag] = fminbnd(@minimizaFocus,-1,1,options); % minimiza en un rango de +/-1 milímetros, el rango es pequeño porque con rangos de +/-5mm se quedaba en mínimos locales ¿?¿
                    % habría que recuperar la posición del motor si no encuentra el foco
                    if (exitFlag == 1)  %% focus found
                        focusString = ['Focus index before=' num2str(-focus_antes) ' after=' num2str(-focusVal) ' relative defocus z=' num2str(z)  ];
                        displaysave(focusString);     % focus history is saved in experiment.txt
                        defocus_mm = z; % the new position is the defocusing
                    else
                        % z es la última posición absoluta utilizada por el minimizador
                        fprintf(xyz,[' 0 0 ' num2str(-z,'%4.4f') ' r \r']); % con esto vuelve al cero original
                        displaysave(['Warning! focus (image) not found for point ',num2str(point)]);
                        defocus_mm = deltaZ; % if not autofocus, the displacement is the specified deltaZ
                    end
                end %is sound_stop
            end
            
            if (AUTOFOCUS && strcmp(AUTOFOCUS_METHOD,'laserspot')==1)  % get the reference
                act_pixel = GetXYfocus(theta,vid,pixel_shotX,pixel_shotY);
                
                if (act_pixel == [0 0])
                    displaysave(['Warning! focus not found for point ',num2str(point)]);
                    defocus_mm = deltaZ; % if not autofocus, the displacement is the specified deltaZ
                else
                    defocus_mm = (act_pixel(1) - ref_pixel(1)) / calibrated_pixels_mm_Z;
                end
                %disp('voy a mover en z');pause;
                
                motor_string = [' 0 0 ' num2str(defocus_mm,'%3.2f') ' r \r'];
                fprintf(xyz,motor_string);
                displaysave(motor_string);
                focusString = ['Focus deltaX=' int2str(act_pixel(1)-ref_pixel(1)) ' ' 'deltaY=' int2str(act_pixel(2)-ref_pixel(2)) ' defocus=' num2str(defocus_mm) ];
                displaysave(focusString);     % focus history is save in experiment.txt
            end
            
            if (AUTOFOCUS==0 && hacer_autofocus==1) % performs a requested autofocus
                displaysave('Haciendo autofoco...');
                set(autofocus_button,'Enable','off')
                options = optimset('Display','notify','TolX',0.1,'MaxFunEvals',20);
                [z focusVal exitFlag] = fminbnd(@minimizaFocus,-2,2,options); % minimiza en un rango de +/-2 milímetros, el rango es pequeño porque con rangos de +/-5mm se quedaba en mínimos locales ¿?¿
                focusString = ['defocus Z=' num2str(z) ];
                displaysave(focusString);
                set(autofocus_button,'Enable','on');
                sound(sound3);
                hacer_autofocus=0;
            end

            if (~ ONLY_RECORDING_PATH) % if only_recording_path, there is no lambdas vectors and procesaMAT crashes
                %disp('doing 1 sec pause(1)... should it be done?');
                %pause(1); % new method (v27) return the LAST capture frame, it could be captured WHILE the motors are moving, removed v64 to speed things up, images could be blury but who cares?
                frame = captureFrameGuppy; % image BEFORE the experiment
                %figure(h_focus);
                set(0,'CurrentFigure',h_focus);
                subplot(5,1,1);
                if ~NO_SHOW_IMAGE
                    set(h_ims_afocus,'CData',frame); % new! v28 update image without a new imshow, promise speed things up
                    %{
                    %v65 removed, it slows everything and not needed, this
                    redraw was the origin of the lag problem in preview
                    hold on; target(pixel_shotX,pixel_shotY); % draw a cross around
                    if (ESCALA==1) % draw 1mm ticks starting at the laser spot
                        if calibrated_pixels_mm_XY_X > 300   % new v58
                            mm1o01=calibrated_pixels_mm_XY_X*0.1; % line every 100um
                            numberOflines=7;
                        else
                            mm1o01=calibrated_pixels_mm_XY_X; % line every 1mm
                            numberOflines=3;
                        end
                        hold on;
                        for horticks=1:numberOflines % cortos hacia la izquierda
                            line([pixel_shotX-(horticks-1)*mm1o01  pixel_shotX-(horticks-1)*mm1o01],[0 300],'Color','w');
                            line([pixel_shotX-(horticks-1)*mm1o01+1  pixel_shotX-(horticks-1)*mm1o01+1],[0 300],'Color','k');
                        end
                        for horticks=1:numberOflines  % largos hacia la derecha
                            line([pixel_shotX+(horticks)*mm1o01  pixel_shotX+(horticks)*mm1o01],[0 600],'Color','w');
                            line([pixel_shotX+(horticks)*mm1o01+1  pixel_shotX+(horticks)*mm1o01+1],[0 600],'Color','k');
                        end
                    end
                    %}
                end % if NO_SHOW_IMAGE
                %[umem,smem] = memory;  displaysave([ 'used memory (GB):   ' num2str(umem.MemUsedMATLAB/1e9)]);
                %imshow(frame,'InitialMagnification', 67); hold on; target(pixel_shotX,pixel_shotY); % draw a cross around
            end   %      v39 se actualiza la foto siempre
            
        end
        
        if (SOUND && ~skipThisPoint)
            record(recorder);
        end
        
        if (ONLY_RECORDING_PATH==0 && ~skipThisPoint)
            
            if exist('extraCommands')
                if not(strcmp(extraCommands,''))
                    try
                        eval(extraCommands);
                    catch
                        disp('extraCommands gave error!!!');
                    end
                end
            end
            
            proceedAfterCheckSurface = true;
            %new v55 checkSurface
            if CHECK_SURFACE > 0 && CAMERA
                pause(0.3); % TO DO: wait for no movement, I don't know if that command exists in corvus
                frame = captureFrameGuppy;
                matrixSpot = frame(pixel_shotY-2:pixel_shotY+2,pixel_shotX-2:pixel_shotX+2);
                brightnessSpot = mean(mean(matrixSpot));
                if brightnessSpot < checkSurface_black
                    proceedAfterCheckSurface = false; % do not do anything
                else
                    checkSurfaceAlreadyPositive = true;
                    checkSurface_matrix(param1Index,param2Index) = 1; % measured
                end
            end
            
            if CHECK_SURFACE <2 && checkSurfaceAlreadyPositive && ~proceedAfterCheckSurface % if there was at least one measurement, and now the surface disapear, lets move onto next param2, ONLY for check_surface==1
                displaysave('checkSurface negative: breaking to next param2');
                break;
            end
            if proceedAfterCheckSurface
                if CHECK_SURFACE > 0
                    displaysave(['Measuring with brightness = ' num2str(brightnessSpot) ]);
                end
                
                Nwindows = size(MW_centralwl,2);
                if (MULTI_WINDOW==0)
                    Nwindows=1;
                end
                
                for nWindow =1:Nwindows % se ejecuta una sola vez en el modo no MULTI_WINDOW
                    
                    if (MULTI_WINDOW && USE_LOCAL_PIMAX)
                        displaysave([' moviendo espectrómetro a ventana ' num2str(MW_centralwl(nWindow)) ]);
                        % objSpe.GetParam('SPT_CUR_POSITION') devuelve la lambda central
                        objSpes.Current.SetParam('SPT_NEW_POSITION',MW_centralwl(nWindow));
                        objSpes.Current.Move();
                        mw_filename = strcat (filename,'_W',sprintf('%02i',nWindow));
                        objExp.SetParam('EXP_DATFILENAME',strcat(folderPath,experimentPath,datePath,mw_filename,'.spe'));
                    end
                    
                    % if using dual pulse laser (software-controlled shutter and qswitch) but also the stm32 pulse generator, we prefer to use the Qnnn feature of the pulse generator
                    % we open the shutter and qswitches now
                    if (USE_DUAL_PULSE_LASER && USE_PULSE_GENERATOR)
                        fwrite(lotis,['SHUT1ON' char(13)]);
                        fwrite(lotis,['QSW1ON' char(13)]);
                        if (DualPulse) % also channel 2
                            fwrite(lotis,['SHUT2ON' char(13)]);
                            fwrite(lotis,['QSW2ON' char(13)]);
                        end
                        pause(0.05);
                        if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
                    end
                    if (USE_SHUTTER_CTRL)
                        %fprintf(shutter,'ens?');pause(0.2);if(shutter.BytesAvailable>0);response=fread(shutter,shutter.BytesAvailable);end
                        %if (numel(response) >= 6 && response(6)=='1')
                        fprintf(shutter,'ens');pause(0.2);if(shutter.BytesAvailable>0);response=fread(shutter,shutter.BytesAvailable);end  % open the shutter
                        %end
                    end
                    
                    % if DEFOCUS we move the sample in Z before and after
                    if ( DEFOCUS ~= 0)
                        displaysave(['DefFocus Z=' num2str(DEFOCUS) ]);
                        moveString = sprintf(' 0 0 %0.3f r \r',DEFOCUS);  fprintf(xyz,moveString); %move axis relative
                    end
                    % if SCANLINE we need to start moving before shooting
                    if ( SCANLINE > 0 )
                        dxSL =  scanLineLength/2 * cos(pi/180*scanLineAngle);
                        dySL =  scanLineLength/2 * sin(pi/180*scanLineAngle);
                        scanLineVel = scanLineLength / Npulses * frequencyPulses;
                        displaysave(['ScanLine first movement XY ' num2str(dxSL) ' ' num2str(dySL) ' for length(mm)=' num2str(scanLineLength,'%2.2f') ]);
                        if SCANLINE <3 % new v53 in mode 3, there is no an initial movement halfway back
                            moveString = sprintf(' %0.3f %0.3f %0.3f r \r',dySL, dxSL, 0); fprintf(xyz,moveString); %move axis relative
                        end
                        % speed change and start moving
                        moveString = sprintf(' %0.3f sv \r',scanLineVel);  fprintf(xyz,moveString); %move axis relative
                        moveString = sprintf(' %0.3f %0.3f %0.3f r \r',-2*dySL, -2*dxSL, 0); fprintf(xyz,moveString); %move axis relative
                        % for next iteration, apply scanLineShrink
                        scanLineLength = scanLineLength * scanLineShrink;
                    end
                    
                    if (USE_ARDUINO_QSWITCH)  fprintf(qswitch,CMD_QSWITCH);    end
                    if (USE_REMOTE_PIMAX==1)
                        try
                            sock = msconnect(remote_pimax_ip,3000);
                            mssend(sock,'measure',2); % 2=timeout
                        catch
                            disp('Error TCP/IP accediendo a REMOTE_PIMAX - measure');
                        end
                    end
                    %time_pulsegen = tic();   xxxxxxx
                    if (USE_PULSE_GENERATOR==1) fprintf(pulsegen,CMD_QSWITCH); end % activates qswitch pulses with NpulsesClean no-qsw pulses
                    
                    %%%%%%%%%%% MEASUREMENT!!  not-autocollect
                    if (USE_LOCAL_PIMAX)
                        objExp.Start(objDoc);
                    end
                    
                    startCapturing = tic();
                    alreadyQswitchAct = 0;
                    % measurement will start shortly, let's check how much time has elapsed from last motor movement
                    displaysave(['  ' datestr(now) ' > time from movement to measure : ' num2str(toc(justAbsMovDone)) ]);
                    for n_spectra=1:Npulses % in theory we must collect Npulses spectra
                        
                        %Qswitch should be activated after NpulsesClean+2 pulses, only if
                        %we are not using the custom stm32 pulse generator
                        if( alreadyQswitchAct==0 && ~USE_PULSE_GENERATOR)
                            if( USE_DUAL_PULSE_LASER  && (toc(startCapturing) > (1 / frequencyPulses * (NpulsesClean+2))))
                                fwrite(lotis,['SHUT1ON' char(13)]);
                                fwrite(lotis,['QSW1ON' char(13)]);
                                if (DualPulse) % also channel 2
                                    fwrite(lotis,['SHUT2ON' char(13)]);
                                    fwrite(lotis,['QSW2ON' char(13)]);
                                end
                                pause(0.05);
                                if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
                                alreadyQswitchAct =1;
                            end
                        end
                        %time_pulsegen_toc = toc(time_pulsegen);
                        %disp(['time1: ' num2str(time_pulsegen_toc)]);
                        
                        if(USE_OCEANOPTICS) % collect Npulse spectra from HR2000 in hardware-external-trigger
                            hr2000_spectrum = wrapper.getSpectrum(0);
                            nast=floor((sum(hr2000_spectrum)-2200000)/20000);
                            nast=max(hr2000_spectrum)/100-10;
                            ast=strcat(num2str(nast),' >');
                            for n=1:nast ast=strcat(ast,'*'); end
                            disp(ast);
                            if (nehr2==10)
                                nast=floor((sum(hr2000_spectrum)-2200000)/20000);
                                displaysave(['Intensidad total HR2000=' num2str(nast)]);
                            end
                            hr2000_spectra(:,n_spectra,param1Index,param2Index) = hr2000_spectrum;
                        end % use_hr2000
                        if(USE_AVANTES) % collect Npulse spectra from Avantes in hardware-external-trigger
                            spectrometer('measure',h1,1); spectrometer('measure',h2,1);  spectrometer('measure',h3,1); spectrometer('measure',h4,1);  %
                            spectrometer('measure',h5,1); spectrometer('measure',h6,1);  spectrometer('measure',h7,1); spectrometer('measure',h8,1);  %
                            waitingForReady = tic();
                            ready=false;
                            while (ready==false)
                                if ready==false
                                    ready=spectrometer('ready',h8); % h5 es el master, aunque con el cable pin6 da igual
                                end
                                pause(0.001); %seconds !!
                                if toc(waitingForReady) > Npulses * 0.1 * 1.5 % a 50% more time than expected means something bad, really bad
                                    panicAvantes = 1;
                                    finish=0; % spectrometer failed on 23feb2023, changed to "=0" to check if it can recover from it
                                    ready=1; % not really true, but let see...
                                    counterAvantesErrors = counterAvantesErrors +1;
                                end
                            end
                            pause(0.01); % just in case not all channels are ready
                            if ~panicAvantes
                                myData1 = spectrometer('getdata',h1);myData2 = spectrometer('getdata',h2);myData3 = spectrometer('getdata',h3);myData4 = spectrometer('getdata',h4);myData5 = spectrometer('getdata',h5);myData6 = spectrometer('getdata',h6);myData7 = spectrometer('getdata',h7);myData8 = spectrometer('getdata',h8);
                                if STORE_SPECTRA8RAW %new v58
                                    spectra8raw(:,n_spectra,param1Index,param2Index,1) = myData1;spectra8raw(:,n_spectra,param1Index,param2Index,2) = myData2;spectra8raw(:,n_spectra,param1Index,param2Index,3) = myData3;spectra8raw(:,n_spectra,param1Index,param2Index,4) = myData4;spectra8raw(:,n_spectra,param1Index,param2Index,5) = myData5;spectra8raw(:,n_spectra,param1Index,param2Index,6) = myData6;spectra8raw(:,n_spectra,param1Index,param2Index,7) = myData7;spectra8raw(:,n_spectra,param1Index,param2Index,8) = myData8;
                                end
                                %all_data = cat(1,myData1(1:1921),myData2(1:1867),myData3(1:1786),myData4(1:1755), myData5(17:1943),myData6(1:1859),myData7(1:2048),myData8(1:2048));
                                all_data= cat(1,myData1(1:1921),myData2(1:1867),myData3(1:1682),myData4(1:1694),myData5(1:1930),myData6(1:1859),myData7(1:2048),myData8(1:2048));
                                if (size(all_data,1) ~= size(avantes_lambdas,1))
                                    error('Wavelength calibration has changed, update concatenation');
                                end
                                %avantes_spectra(:,n_spectra,param1Index,param2Index) = all_data;
                                spectra(:,n_spectra,param1Index,param2Index) = all_data;
                            else
                                if counterAvantesErrors > 5 
                                   finish=1;
                                end
                            end % panic
                        end % use_avantes
                        
                        if (USE_LOCAL_PIMAX)
                            % do nothing, it is measuring in the background
                        end % USE_LOCAL_PIMAX
                        
                    end % for npulse
                    
                    % if only USE_LOCAL_PIMAX is used, it will arrive here too soon
                    if (USE_LOCAL_PIMAX)
                        while objExp.GetParam('EXP_RUNNING')
                            % waiting for the experiment to finish
                        end
                    end
                    if (USE_PULSE_GENERATOR==1) fprintf(pulsegen,['q0' char(13)]); end % deactivates qswitch pulses
                    if (USE_DUAL_PULSE_LASER)
                        fwrite(lotis,['QSW1OFF' char(13)]);
                        fwrite(lotis,['SHUT1OFF' char(13)]);
                        if (DualPulse) % also channel 2
                            fwrite(lotis,['QSW2OFF' char(13)]);
                            fwrite(lotis,['SHUT2OFF' char(13)]);
                        end
                        pause(0.05);
                        if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
                    end
                    if (USE_ARDUINO_QSWITCH)  fprintf(qswitch,'0');    end
                    if (USE_SHUTTER_CTRL)
                        fprintf(shutter,'ens');pause(0.2);if(shutter.BytesAvailable>0);response=fread(shutter,shutter.BytesAvailable);end  % close the shutter
                    end
                    
                    if (USE_LOCAL_PIMAX)
                        objDoc.Save;  % save the datafile
                    end
                    
                    if (USE_REMOTE_PIMAX==1)
                        try
                            ltemp = msrecv(sock,2); % lambdas
                            stemp = msrecv(sock,2); % spectra
                            ack = msrecv(sock,2);  % ok
                            lambdas_remote_pimax = ltemp;
                            spectra_remote_pimax(:,:,param1Index,param2Index) = stemp;
                            msclose(sock);
                        catch
                            disp('Error TCP/IP accediendo a REMOTE_PIMAX - espectro');
                        end
                    end
                    
                    % mar2017 trying to give time to winspec to do its things with the temporal file, just in case
                    if (~ONLY_RECORDING_PATH)
                        %disp('doing 1 sec pause(2)... should it be done?');
                        %pause(1);
                    end
                    
                    % if SCANLINE we need to go back to the point at warp speed
                    if ( SCANLINE > 0)
                        fprintf(xyz,' 15 sv \r'); %move axis relative
                        if SCANLINE < 3 % 1 & 2 mode: half length back
                            moveString = sprintf(' %0.3f %0.3f %0.3f r \r',dySL, dxSL, 0); fprintf(xyz,moveString); %move axis relative
                        else % new v53 no initial movement, we need to go back the entire lentgh
                            moveString = sprintf(' %0.3f %0.3f %0.3f r \r',2*dySL, 2*dxSL, 0); fprintf(xyz,moveString); %move axis relative
                        end
                        if (scanLineLength < 0.02) % too small, %v33 se peta el motor por ser numeros pequeños
                            SCANLINE=0;
                        end
                    end
                    if ( DEFOCUS ~= 0) % the same movement but in oposite direction
                        moveString = sprintf(' 0 0 %0.3f r \r',- DEFOCUS);  fprintf(xyz,moveString); %move axis relative
                    end
                    
                    if (first_time_lambdas==1)  % only with PIMAX we need to get the lambdas vector after measurement
                        if (USE_LOCAL_PIMAX)
                            objCal = objDoc.GetCalibration; % get calibration data from the last experiment
                            % create array of lambdas
                            if (MULTI_WINDOW==1)
                                for x = 1:Npixeles
                                    MW_lambdas(x,nWindow) = objCal.Lambda(x);
                                end
                                % en MULTI_WINDOW estas lambdas pueden estar repartidas en var w
                                for x = 1:size(satlines,2)
                                    p=buscaPixel(MW_lambdas(:,nWindow),satlines(x)); % no hay espectro pero da igual porque no se busca el máximo real
                                    if (p>1 && p<1024) % lo ha encontrado
                                        satlines_pix(x) = p;
                                    end
                                end
                            else
                                for x = 1:Npixeles
                                    lambdas(x) = objCal.Lambda(x);
                                end
                            end
                        else
                            lambdas = avantes_lambdas;
                        end % USE_LOCAL_PIMAX
                        
                        for x = 1:size(satlines,2)
                            satlines_pix(x) =buscaPixel(lambdas,satlines(x)); % no hay espectro pero da igual porque no se busca el máximo real
                        end
                    end
                    
                    % read back the spectral information from the just saved file using
                    % readSPE.m
                    if (USE_LOCAL_PIMAX)
                        try % sometimes winspec stop with error "read only tmp file"
                            objDoc.Close; % close the window
                            if (MULTI_WINDOW)
                                tmp_spectra = readSPE(strcat(folderPath,experimentPath,datePath,mw_filename,'.spe'));
                                MW_spectra(:,1:size(tmp_spectra,3),param1Index,param2Index,nWindow) = tmp_spectra(1,:,:);  % zeros(Npixeles,Npulses, Nparam1, Nparam2,size(MW_centralwl,2),'single');
                            else
                                tmp_spectra = readSPE(strcat(folderPath,experimentPath,datePath,filename,'.spe'));
                                spectra(:,1:size(tmp_spectra,3),param1Index,param2Index) = tmp_spectra(1,:,:);
                            end
                        catch
                            displaysave('WinSpec Error with temporal files!!! hay que volver a cargar winSpec');
                            sound(sound2);sound(sound3); % sound alert
                            pause;
                            
                            %% WixX32 objects
                            objCal = actxserver('WinX32.CalibObj');
                            objWin = actxserver('WinX32.DocWindows');
                            objDocs = actxserver('WinX32.DocFiles');    % It seems needed to access the calibration of the first DOC
                            objDoc = actxserver('WinX32.DocFile');    %
                            objExp = actxserver('WinX32.ExpSetup');
                            objSpe = actxserver('WinX32.SpectroObj');
                            objPul = actxserver('WinX32.PITG');
                            objSpes = actxserver('WinX32.SpectroObjMgr');
                            % gain
                            actual_gainIntensifier = gainIntensifier;
                            objExp.SetParam('EXP_INTENSIFER_GAIN',actual_gainIntensifier);
                            % pulser
                            objPul.SetParam('TGC_PULSE_DELAY',delayCapture);
                            objPul.SetParam('TGC_TRIG_FREQUENCY',frequencyPulses);
                            objPul.SetParam('TGC_PULSE_WIDTH',gateWidth);
                            objPul.Process('TGP_SET_GATE_PULSE'); % commits the software parameters into hardware
                        end
                    end
                    
                    % plot the largest spectra
                    %figure(h_focus);
                    set(0,'CurrentFigure',h_focus);
                    subplot(5,1,2); % spectra
                    max_sum = 0;
                    cual=1;
                    for n=1:Npulses
                        energia = sum(spectra(:,n,param1Index,param2Index));
                        if ( energia > max_sum)
                            max_sum = energia;
                            cual = n;
                        end
                    end
                    plot(lambdas,spectra(:,cual,param1Index,param2Index));
                    axis([ 150 950 0 66000]);
                    energia = (sum(spectra(:,cual,param1Index,param2Index)) - 13E6) / 1E5;
                    if (energia > energia_limit)
                        energia_limit = energia_limit * 10;
                    end
                    
                    %%OJO!! the energy bar has been substituted by the remote pimax
                    %%spectrum
                    %figure(h_focus);
                    set(0,'CurrentFigure',h_focus);
                    subplot(5,1,3);  %bar
                    if (USE_REMOTE_PIMAX)
                        plot(ltemp,stemp);
                        axis([ ltemp(1) ltemp(end) 0 66000]);
                    else
                        barh(energia);
                        xlim([0 energia_limit]);
                    end
                    
                    %v55 new figure at subplot(5,1,5) with the lineMonitoring information
                    h_focus_lines = barh(lineMonHeights);
                    set(gca,'XScale','log');
                    xlim([800 65000]); % new v58 raw value to better detect saturation
                    yticklabels(lineMonitoring_name); % absolute energy in exponential format
                    hay_red = false;
                    hay_blue = false;
                    lineMonString='';
                    for i=1:numel(lineMonitoring_lambda)
                        %new v58 search the peaks
                        idx= lambda2pixel(avantes_lambdas,lineMonitoring_lambda(i));
                        [maxval maxidx] = max(spectra(idx-5:idx+5,cual,param1Index,param2Index));
                        newidx = idx-5+maxidx-1;
                        %lineMonitoring_lambda(i) = avantes_lambdas(newidx); no need to change
                        val = spectra(newidx,cual,param1Index,param2Index); % peak height
                        if val > ValSaturation
                            lineMonString = strcat(lineMonString,'  SAT!',lineMonitoring_name(i),':',num2str(val,'%05i ')); % magenta if value is higher than max
                            hay_red=true;
                            displaysave('WARNING! saturation of monitored line!');
                        elseif (val>lineMonitoring_max(i))
                            lineMonString = strcat(lineMonString,'  ^',lineMonitoring_name(i),':',num2str(val,'%05i ')); % magenta if value is higher than max
                        elseif (val < 0.5*lineMonitoring_min(i))
                            lineMonString = strcat(lineMonString,'  .',lineMonitoring_name(i),':',num2str(val,'%05i ')); % red if it is  lower than min
                            hay_blue=true;
                        else
                            lineMonString = strcat(lineMonString,'  =',lineMonitoring_name(i),':',num2str(val,'%05i ')); % just fine
                        end
                        lineMonHeights(i) = val;
                    end
                    h_focus_lines.YData = lineMonHeights;
                    if hay_red
                        set(h_focus_lines,'FaceColor','r'); % at least one line is saturated
                    elseif hay_blue
                        set(h_focus_lines,'FaceColor','m'); % at least one line below 50% of minimum reference value
                    else
                        set(h_focus_lines,'FaceColor','b');
                    end
                    refreshdata(h_focus_lines,'caller');
                    
                    %% new v57 temperature recording of avantes spectrometer internal thermistor
                    if (USE_AVANTES)
                        ther = spectrometer('getanalogin',h5,6);
                        spectrometerTemperature(param1Index,param2Index) = 118.69- 70.361*ther+21.02*ther*ther-3.6443*ther*ther*ther+0.1993*ther*ther*ther*ther;
                    end
                    
                    %% processing after a complete point/experiment
                    
                    if (SOUND) % save audio data to a file
                        stop(recorder);
                        %function audiowrite is not available in R2010a
                        %audiowrite(strcat(FolderPath,experimentPath,datePath,filename,'.wav'),getaudiodata(recorder),44100);
                        tmp_audio = getaudiodata(recorder);
                        %save(strcat(folderPath,experimentPath,datePath,filename,'_audio','.mat'),'tmp_audio');
                        %wavwrite(tmp_audio,44100,strcat(folderPath,experimentPath,datePath,filename,'_audio.wav'));
                        audiowrite(strcat(folderPath,experimentPath,datePath,filename,'_audio.wav'),tmp_audio,44100,'BitsPerSample',24);
                        
                        if (SOUND_STOP ==1 ) % monitorizado del sonido
                            %proceso el archivo de sonido
                            bhi = fir1(34,0.1,'high',chebwin(35,30)); % copiado de un ejemplo, no está claro algún parámetro, la frecuencia es el 10%
                            ondaf = filter(bhi,1,tmp_audio);
                            r=rms(ondaf,50,30,0); % la función rms está descargada y es "windowed"
                            [p,idx]=findpeaks(r,'MINPEAKHEIGHT',max(r)/5);
                            % ahora hay que calcular la energía de cada pulso
                            integrated_sound = zeros (size(p,2),1);
                            distancia = round((idx(2)-idx(1))/10);
                            % ha dado un error porque distancia = 1083 y se sale del borde
                            for j=1:(length(p))  % OJO!! que puse un -1 y no sé por qué, integrated_sound tiene así un cero al final
                                integrated_sound(j) = sum(r(idx(j)-distancia:idx(j)+distancia));
                            end
                            energia_sonido_promedio = mean(integrated_sound);
                            %std_sonido = movingstd(integrated_sound,10); k tiene que ser
                            %length(p)/2-1
                            std_sonido = std(integrated_sound);
                            std_sonido_promedio = mean(std_sonido);
                            clear tmp_audio; % just to save memory
                            clear ondaf;
                            clear r;
                            displaysave(['Pulsos de sonido: ' num2str(length(p)) ' Energía: ' num2str(energia_sonido_promedio)  ' STD acústica: ' num2str(std_sonido_promedio)]);
                            
                            if (param1Index==1) %% no hago nada hasta tener dos medidas
                                std_sonido_ref1 = std_sonido_promedio;
                            elseif (param1Index==2) %% no hago nada hasta tener dos medidas
                                std_sonido_ref2 = std_sonido_promedio;
                            elseif (param1Index==3) %% hay dos medidas
                                std_sonido_referencia =  (std_sonido_ref1 + std_sonido_ref2)/2;
                            else
                                if (std_sonido_promedio > STD_SONIDO_STOP*std_sonido_referencia) % un xx% más, suele ser muy estable al principio
                                    sound(sound2);sound(sound3); % sound alert
                                    disp('Posible perforación por el sonido... haciendo skipping=1');
                                    skipping=1;
                                end
                            end % if param1Index==1
                        end
                    end % if sound
                end % multi_window . se ejecuta una sola vez si MULTI_WINDOW=0
                
                first_time_lambdas=0;
                pointCounterEffective = pointCounterEffective +1; % new v14/v64
                
            else % proceedAfterCheckSurface
                displaysave(['Skipping checkSurface with brightness = ' num2str(brightnessSpot) ' at param1=' num2str(param1) ',param2=' num2str(param2)]);
            end % if checkSurface
            
            %% AUTOGAIN - NOAUTOCOLLECT. if the mean intensity of the monitored lines outside a +/-50% range , a new gain is calculated. if the new gain is > MAX_GAIN or <1, does not change.
            if (AUTOGAIN==1) % old algorithm, look into peak saturation
                gain_intensity = mean (tmp_spectra(1,satlines_pix,:),3);
                mean_intensity = mean(gain_intensity);
                pixeles_saturados = find(tmp_spectra);
                
                if (num_broken > Npulses/20) % if 20% of spectra are broken
                    actual_gainIntensifier = floor (actual_gainIntensifier/2);
                    if (actual_gainIntensifier > 0 && actual_gainIntensifier < MAX_GAIN)
                        objExp.SetParam('EXP_INTENSIFER_GAIN',actual_gainIntensifier);
                        toprint= [' >Gain changed due to ' num2str(num_broken) ' broken spectra to ', num2str(actual_gainIntensifier), ' after point ', num2str(point)];
                        displaysave(toprint);
                    end
                    
                else  if (mean_intensity > 1.5*MEAN_INTENSITY || mean_intensity < 0.5*MEAN_INTENSITY)
                        actual_gainIntensifier = floor (MEAN_INTENSITY / mean_intensity * actual_gainIntensifier);
                        if (actual_gainIntensifier > 0 && actual_gainIntensifier < MAX_GAIN)
                            objExp.SetParam('EXP_INTENSIFER_GAIN',actual_gainIntensifier);
                            toprint= [' >Gain changed to ', num2str(actual_gainIntensifier), ' after point ', num2str(point)];
                            displaysave(toprint);
                        end
                    end
                    autoGains(param1Index,param2Index) = actual_gainIntensifier;
                end
            elseif (AUTOGAIN==2)  % try to maintain the mean global intensity of the first measurement
                if (MULTI_WINDOW)
                    spectra = MW_spectra(:,:,:,:,1); % se usa la primera ventana
                end
                first_intensity = sum(sum(spectra(:,NpulsesClean:end,1,1),2),1);
                first_IG_ratio = first_intensity * gainIntensifier; % original intensity x gain ratio, to be maintained
                actual_intensity = sum(sum(spectra(:,NpulsesClean:end,param1Index,param2Index)));
                actual_gainIntensifier = round( first_IG_ratio / actual_intensity * 1.4) ; % adhoc factor because gain factor does not traslate to amplitude increment
                if (actual_gainIntensifier-previous_gain > 3) % capped
                    actual_gainIntensifier = previous_gain+3;
                end
                if (actual_gainIntensifier-previous_gain < 0) % no se permita que baje
                    actual_gainIntensifier = previous_gain;
                end
                if (actual_gainIntensifier > 80)
                    actual_gainIntensifier=80; % jut to not petarlo
                end
                objExp.SetParam('EXP_INTENSIFER_GAIN',actual_gainIntensifier);
                toprint= [' >AutoGain changed to ' num2str(actual_gainIntensifier) ' after point ' num2str(point) ' with integrated intensity ' num2str(actual_intensity) ' first intensity=' num2str(first_intensity) ' first_IG_ratio=' num2str(first_IG_ratio) ];
                displaysave(toprint);
                autoGains(param1Index,param2Index) = actual_gainIntensifier;
                previous_gain = actual_gainIntensifier;
            end
            
            
            
            %% Shows statistics of the experiment
            if (MULTI_WINDOW)
                spectra = MW_spectra(:,:,:,:,1); % se usa la primera ventana
                lambdas = MW_lambdas(:,1);
            end
            num_withoutpeaks=0;
            num_saturated=0;
            num_valid_spectra=0;
            num_broken=0;
            is_withoutpeaks = 0;
            is_saturated =0;
            is_broken=0;
            for i=1:Npulses
                % primero ver los picos en los que detectar saturación
                if (size(satlines,2) == 0 || min(min(satlines_pix,1))==0)  % there is no specific line to watch for or outside the window
                    peak_max = max(spectra(:,i,param1Index,param2Index));
                else % only the satlines lines are compared to ValSaturation
                    peak_max = 0;
                    for x=1:size(satlines,2)
                        if (spectra(satlines_pix(x),i,param1Index,param2Index) > peak_max)
                            peak_max = spectra(satlines_pix(x),i,param1Index,param2Index);
                        end
                    end
                end
                peak_mean = mean(spectra(:,i,param1Index,param2Index));
                spectra_is_valid(i,param1Index,param2Index)=1; % for the moment...
                if (size(find(spectra(:,i,param1Index,param2Index)<10),1) > 0) % if there is any intensity value lower than 10, the spectrum is "broken" (excess of light)
                    num_broken = num_broken+1;
                    is_broken=0;
                end
                
                if (peak_max/peak_mean < 0)   % OJO!!! está anulado porque puedo monitorizar lìneas de magnesio muy flojas
                    num_withoutpeaks = num_withoutpeaks +1;
                    is_withoutpeaks=1;
                end
                if (peak_max > ValSaturation)
                    num_saturated = num_saturated +1;
                    is_saturated=1;
                end
                if (is_broken || is_saturated || is_withoutpeaks)
                    spectra_is_valid(i,param1Index,param2Index)=0; % what a pity
                else
                    num_valid_spectra = num_valid_spectra +1;
                end
            end
            
            num_valid_spectra_matrix(param1Index,param2Index) = num_valid_spectra;
            intensidad = sum(sum(spectra(:,:,param1Index,param2Index),2),1);
            toprint=['  >', num2str(Npulses),' spectra,  ',num2str(num_valid_spectra), ' valids,  ',num2str(num_saturated), ' saturated,  ',num2str(num_withoutpeaks),' without peaks, ',num2str(num_broken),' broken.' ' Total intensity: ' num2str(intensidad, '%1.2e')];
            displaysave(toprint);
            
            % print heights of lines under monitoring
            lineMonString = strings;
            tmp_mean_spectra = mean(spectra(:,NpulsesClean+2:end,param1Index,param2Index),2); % won't work for only one or two shots
            for i=1:numel(lineMonitoring_lambda)
                val = tmp_mean_spectra(lambda2pixel(avantes_lambdas,lineMonitoring_lambda(i)));
                if (val>lineMonitoring_max(i))
                    lineMonString = strcat(lineMonString,'  ^',lineMonitoring_name(i),':',num2str(val,'%05i ')); % magenta if value is higher than max
                elseif (val < lineMonitoring_min(i))
                    lineMonString = strcat(lineMonString,'  .',lineMonitoring_name(i),':',num2str(val,'%05i ')); % red if it is  lower than min
                else
                    lineMonString = strcat(lineMonString,'  =',lineMonitoring_name(i),':',num2str(val,'%05i ')); % just fine
                end
                lineMonHeights(i) = val/((lineMonitoring_max(i)+lineMonitoring_min(i))/2)*100; % 100% = (max+min)/2
            end
            displaysave(lineMonString);
            
            
            
            if (AUTOGAIN && num_saturated > Npulses/10) % if more than 10% of spectra are saturated
                actual_gainIntensifier = floor(actual_gainIntensifier * 0.9);
                objExp.SetParam('EXP_INTENSIFER_GAIN',actual_gainIntensifier);
                toprint= ['   >Gain reduced to ', num2str(actual_gainIntensifier), ' after point ', num2str(point)];
                displaysave(toprint);
            end
            
        end % if only_recording_path
        
        
        if (RSD_PROCESSING==1 && ONLY_RECORDING_PATH==0)
            if (MULTI_WINDOW)
                spectra = MW_spectra(:,:,:,:,1); % se usa la primera ventana
                lambdas = MW_lambdas(:,1);
            end
            
            %RSD calc at one single pixel
            onlyvalids = find(spectra_is_valid(:,param1Index,param2Index)==1); % the indexes of valid spectra
            tmp_spectra_mean_valids = nanmean(spectra(:,onlyvalids,param1Index,param2Index),2);
            
            if (pixel_captured_from_first_spectrum==0)
                pixelCa = buscaPixel(lambdas,LAMBDA_CALCIO,tmp_spectra_mean_valids,0); % busca el pixel con flag=0, se supone que hay picos porque están seleccionados los válidos
                pixelMg = buscaPixel(lambdas,LAMBDA_MAGNESIO,tmp_spectra_mean_valids,0);
                pixelBackGroundCentral_Mg = buscaPixel(lambdas,(background4lambdas_mg(1)+background4lambdas_mg(2)) / 2,tmp_spectra_mean_valids,2); % punto central del background izquierdo
                pixelBackGroundCentral_Ca = buscaPixel(lambdas,(background4lambdas_ca(1)+background4lambdas_ca(2)) / 2,tmp_spectra_mean_valids,2);
                pixel_captured_from_first_spectrum=1;
            end
            
            all_intensities(:,param1Index,param2Index) = sum(spectra(:,:,param1Index,param2Index),1); % the total intensities of ALL the Npulses spectra
            tmp_spectra_rsd = spectra(pixelCa,onlyvalids,param1Index,param2Index);
            tmp_spectra_rsdratio = spectra(pixelMg,onlyvalids,param1Index,param2Index) ./ spectra(pixelCa,onlyvalids,param1Index,param2Index);
            
            sRSD = std(tmp_spectra_rsd); %standard deviation of a single point and a single pixel across the Npulses
            sInt = std(all_intensities(:,param1Index,param2Index)); % std of the integrated intensisties
            mInt = mean(all_intensities(:,param1Index,param2Index)); % mean  of the integrated intensisties
            sRSDratio = std(tmp_spectra_rsdratio);
            mRSD = mean(tmp_spectra_rsd);%mean to get the RELATIVE std
            mRSDratio = mean(tmp_spectra_rsdratio);
            ratio_MgCa_toplot(param1Index,param2Index) = mRSDratio;
            
            % calculo de SNR
            tmp_spectra_snr = spectra(:,onlyvalids,param1Index,param2Index);
            
            intensidad_Mg = tmp_spectra_snr(pixelMg,:); % intensidades del pixel del pico de magnesio y calcio
            intensidad_Ca = tmp_spectra_snr(pixelCa,:);
            intensidad_backg_Mg = tmp_spectra_snr(pixelBackGroundCentral_Mg,:); % todas las intensidades de los espectros, en el pixel del background, para hacer luego el std
            intensidad_backg_Ca = tmp_spectra_snr(pixelBackGroundCentral_Ca,:);
            SNR_Mg(param1Index,param2Index) = nanmean(intensidad_Mg) / nanstd(intensidad_backg_Mg);
            SNR_Ca(param1Index,param2Index) = nanmean(intensidad_Ca) / nanstd(intensidad_backg_Ca);
            
            
            %toprint=[' Peak at ',sprintf(' %.2f',lambdas(pixelCa)),'nm  mean=',sprintf(' %.2f',mRSD),' rsd Ca line=',sprintf(' %2.2f',100*sRSD/mRSD),'%%',' Mg/Ca ratio mean=',sprintf(' %.2f',mRSDratio),' Mg/Ca ratio rsd=',sprintf(' %2.2f',100*sRSDratio/mRSDratio),'%%',' all_intensities rsd=',sprintf(' %2.2f',100*sInt/mInt),'%%'];
            %deleted v55 toprint=[' Mg/Ca ratio mean=',sprintf(' %.2f',mRSDratio),' Mg/Ca ratio rsd=',sprintf(' %2.2f',100*sRSDratio/mRSDratio),'%%',' all_intensities rsd=',sprintf(' %2.2f',100*sInt/mInt),'%%',' Peak at ',sprintf(' %.2f',lambdas(pixelCa)),'nm  mean=',sprintf(' %.2f',mRSD),' rsd Ca line=',sprintf(' %2.2f',100*sRSD/mRSD)];
            %displaysave(toprint);
            %toprint=['SNR_Ca=', sprintf('%3.1f',SNR_Ca(param1Index,param2Index)),' SNR_Mg=',sprintf('%3.1f',SNR_Mg(param1Index,param2Index))];
            %displaysave(toprint);
            
            Rsd_1pixel(param1Index,param2Index) = sRSD/ mRSD;  Rsd_1pixel_mean(param1Index,param2Index)= mRSD;
            Rsd_ratio(param1Index,param2Index) = sRSDratio/ mRSDratio;  Rsd_ratio_mean(param1Index,param2Index)= mRSDratio;
            
            % pinta la secuencia de ratios hasta el momento
            %figure(h_focus);
            set(0,'CurrentFigure',h_focus);
            subplot(5,1,4);
            ylim('auto');
            eje = param1From:param1Step:param1To;
            plot(eje,ratio_MgCa_toplot(:,param2Index));
            %ylim([0 2]);
            
        end
        
        %% camera image after experiment
        if (CAMERA&& ~skipThisPoint)
            if (ONLY_1PHOTO==0 || (ONLY_1PHOTO==1 && param1Index==1 && Nparam2>1)) % only takes photo for the first shoots of every param2
                if (~ ONLY_RECORDING_PATH) % if only_recording_path, no photos to speed up things
                    pause(0.6); % new method (v27) return the LAST capture frame, it could be captured WHILE the motors are moving, v64 0.6sec seems a safe delay even for low light
                    try
                        frame = captureFrameGuppy;
                        imwrite(frame,strcat(folderPath,experimentPath,datePath,filename,'.png')); %save current frame with the same name as the experiment
                    catch
                    end
                end
            end
        end
        
        % imagen de la webcam que debería apuntar al medidor de energía de
        % pulsos
        if (TAKE_WEBCAM_PHOTO && ~skipThisPoint)
            if (~ ONLY_RECORDING_PATH) % if only_recording_path, no photos to speed up things
                frame = captureFrameWebcam;
                try
                    imwrite(frame,strcat(folderPath,experimentPath,datePath,filename,'_webcam.png')); %save current frame with the same name as the experiment
                catch
                end
            end
        end
        
        % aqui estaba el modo automático
        % transverse movement is 2D scanning is activated (param1Name = deltaT)
        if ( strcmp(param1Name,'deltaT') == 1)
            if (param1Index == Nparam1) % last movement
                dx2D =  L2D/2 * cos(pi/180*DeltaT_ANGLE);
                dy2D =  L2D/2 * sin(pi/180*DeltaT_ANGLE);
                displaysave(['deltaT last mov XY ' num2str(dx2D) ' ' num2str(dy2D)]);
            else
                dx2D = - L2D/(Npoints2D-1) * cos(pi/180*DeltaT_ANGLE);
                dy2D = - L2D/(Npoints2D-1) * sin(pi/180*DeltaT_ANGLE);
                displaysave(['deltaT mov XY ' num2str(param1Index) ' ' num2str(dx2D) ' ' num2str(dy2D)]);
            end
            moveString = sprintf(' %0.3f %0.3f %0.3f r \r',dy2D, dx2D, 0); %
            fprintf(xyz,moveString); %move axis relative
        elseif ( strcmp(param1Name,'pulsos') == 1) % en este modo se sigue capturando espectros sin moverse
            
        else
            %"normal" movement (automatic and manual) moved here
            % new v49 "posX" lop names mean absolute coordinates
            %%start code to move motors, now it is duplicated to avoid nested functions
            justAbsMovDone=tic(); % to check if motors are effectively stopped, init here
            if (MOTORS && ~MANUAL_TRACKING && ~skipThisPoint)
                %first, if SCANLINE=2, we need user input
                if(USE_RECORDED_PATH)
                    deltaX=-old_recorded_path(point,2);deltaY=-old_recorded_path(point,1);deltaZ=old_recorded_path(point,3);
                end
                if absoluteMovement == 0
                    moveString = sprintf(' %0.3f %0.3f %0.3f r \r',deltaY, deltaX, deltaZ); % EY!! OJO!! coordiantes changed from physical to image
                    fprintf(xyz,moveString); %move axis relative
                    justAbsMovDone=tic(); % to check if motors are effectively stopped
                    statusXYZ=1; % v64 code to wait for motors to stop, already tested
                    while statusXYZ==1        %is moving,
                        pause(0.05);
                        if(xyz.BytesAvailable>0);fread(xyz,xyz.BytesAvailable);end
                        fprintf(xyz,' status \r');
                        statusXYZ = fscanf(xyz,'%f'); % status word
                    end
                end
                actual_deltaP = deltaP;
                if (deltaP ~= 0.0) % angle motor
                    if (rem(point,everyNpoint)==0) % this is the nth point with possible different displacement
                        actual_deltaP = deltaP + everyIncDec;
                    end
                    if (USE_RECORDED_PATH)
                        actual_deltaP = old_recorded_path(point,5);
                    end
                    total_distance = total_distance + actual_deltaP + sqrt(deltaX^2+deltaY^2);
                    actual_diameter = DIAMETER - param1 / param1To * (DIAMETER - DIAMETER_END );
                    theta_steps = round( actual_deltaP / (3.14159*actual_diameter / 2400));
                    for anglestep=1:theta_steps
                        if (deltaP > 0)
                            %fwrite(theta,'w');pause(0.1); total_ws_steps = total_ws_steps + 1;
                        else
                            %fwrite(theta,'s');pause(0.1); total_ws_steps = total_ws_steps - 1 ;
                        end
                    end
                end
                recorded_path(point,1) = deltaX;recorded_path(point,2) = deltaY;recorded_path(point,3) = deltaZ;recorded_path(point,4) = defocus_mm;recorded_path(point,5) = actual_deltaP;
            end
            
            %% start in manual_tracking mode
            if (MOTORS && MANUAL_TRACKING)
                if(mantrack_next_points == 0)  % there are no automatic points to follow the path
                    if (~USE_RECORDED_PATH) % not using using_recorded_path
                        %sound(sound2);sound(sound3); % sound alert                       
                        if(finish) % added v33 : finish before asking for click
                            displaysave('Breaking from mantracking-preview-window');
                            break; % OJO!! cambiado 20enero2015 para ver si así termina bien
                        end
                        %twitter
                        if (TWITTER==2)
                            try
                                tw.updateStatus(['Autolibs: requiere nueva dirección manual_tracking ' datestr(now) ]);
                            catch
                                TWITTER=0;
                            end
                        end
                        if haveImageStitched  % new v62 uses the stitched image
                            %get current coordinates
                            fprintf(xyz,' p \r1 j');pause(0.05); % get position and enable joystick again
                            current_coords = fscanf(xyz,'%f %f %f'); % this reads the current coordinates x y z
                            %disp(['  Current coords: x=' num2str(current_coords(2)) ',y=' num2str(current_coords(1))]);
                            %absX = current_coords(2);
                            %absY = current_coords(1);
                            %absZ = current_coords(3);
                            figure(h_stitch);
                            h_stitch.Name = 'Select next point with a click - press other buttons or a key to finish';
                            finishStitch=0; % we reuse this to finish manual_tracking
                            %set(esc_button_stitch,'Enable','off');
                            targetmm(current_coords(2), current_coords(1)); % draw a cross around the current coordinates
                            while(1)
                                %[xc, yc, button] = ginput(1); % button returns how was pressed
                                %v08 alternative with drawpoint
                                clear roi_nextpoint; % this way, if ESC is pressed, return empty
                                roi_nextpoint = drawpoint('Color','red');
                                if size(roi_nextpoint.Position,2) < 2 % empty
                                    finish=1;
                                    break;
                                else
                                    xc = roi_nextpoint.Position(1);
                                    yc = roi_nextpoint.Position(2);
                                end
                                distx = (xc-current_coords(2));
                                disty = (yc-current_coords(1));
                                mantrack_dist =  sqrt((distx)^2+(disty)^2); % distance in mm
                                mantrack_arctan = atan(abs(disty)/abs(distx)); % tangent(angle) v50: with non-square pixels, tan should be done in mm
                                if (xc>current_coords(2))
                                    dirX=1;
                                else
                                    dirX=-1;
                                end
                                if (yc>current_coords(1)) %OJO!! coordinates in stitchedImage are not the same as in old mehtod in pixels
                                    dirY=-1;
                                else
                                    dirY=1;
                                end
                                
                                % wee need to calculate how many points there is until next
                                % stop, taking into acount everyNpoints;
                                if (MANUAL_TRACKING==2)
                                    quedan=1;quedadist=mantrack_dist;
                                else
                                    quedan=0;quedadist=0;
                                    while(1)
                                        if (rem(point+quedan,everyNpoint)==0)
                                            quedadist=quedadist+deltaR+everyIncDec;
                                        else
                                            quedadist=quedadist+deltaR;
                                        end
                                        if (quedadist>mantrack_dist)
                                            break; % at this point, 'quedan' is the number of points
                                        else
                                            quedan=quedan+1;
                                        end
                                        
                                    end
                                end % if manual_tracking=1
                                mantrack_next_points = quedan - 1 ;
                                if (mantrack_next_points>-1)
                                    break;
                                end                           
                            end % while(1)

                            
                            
                        else % old "click" window                            
                            h_click=figure('Position',[1 1 scrsz(3) scrsz(4)],'Name','Current image - click for direction of new points'); % create new large window
                            pause(1); % new method (v27) return the LAST capture frame, it could be captured WHILE the motors are moving
                            frame = captureFrameGuppy;
                            imshow(frame,'InitialMagnification', 67); hold on;
                            target(pixel_shotX,pixel_shotY); % draw a cross around
                            if (ESCALA==1) % draw 1mm ticks starting at the laser spot
                                if calibrated_pixels_mm_XY_X > 300   % new v58
                                    mm1o01=calibrated_pixels_mm_XY_X*0.1; % line every 100um
                                    numberOflines=7;
                                else
                                    mm1o01=calibrated_pixels_mm_XY_X; % line every 1mm
                                    numberOflines=3;
                                end
                                hold on;
                                for horticks=1:numberOflines % cortos hacia la izquierda
                                    line([pixel_shotX-(horticks-1)*mm1o01  pixel_shotX-(horticks-1)*mm1o01],[0 300],'Color','w');
                                    line([pixel_shotX-(horticks-1)*mm1o01+1  pixel_shotX-(horticks-1)*mm1o01+1],[0 300],'Color','k');
                                end
                                for horticks=1:numberOflines  % largos hacia la derecha
                                    line([pixel_shotX+(horticks)*mm1o01  pixel_shotX+(horticks)*mm1o01],[0 600],'Color','w');
                                    line([pixel_shotX+(horticks)*mm1o01+1  pixel_shotX+(horticks)*mm1o01+1],[0 600],'Color','k');
                                end
                            end
                            
                            while(1)
                                %para mostrar un video en esto, https://es.mathworks.com/matlabcentral/answers/275268-how-to-use-ginput-with-live-video
                                %v43 MANUAL_TRACKING=2 move motors just to the target point
                                [xc, yc] = ginput(1);
                                distx = (xc-pixel_shotX)/calibrated_pixels_mm_XY_X;
                                disty = (yc-pixel_shotY)/calibrated_pixels_mm_XY_Y;
                                mantrack_dist =  sqrt((distx)^2+(disty)^2); % distance in mm
                                mantrack_arctan = atan(abs(disty)/abs(distx)); % tangent(angle) v50: with non-square pixels, tan should be done in mm
                                %disp(['distancia: ' num2str(mantrack_dist) ' angulo: ' num2str(180/3.1415*mantrack_arctan)])
                                if (xc>pixel_shotX)
                                    dirX=1;
                                else
                                    dirX=-1;
                                end
                                if (yc>pixel_shotY)
                                    dirY=1;
                                else
                                    dirY=-1;
                                end
                                % wee need to calculate how many points there is until next
                                % stop, taking into acount everyNpoints;
                                if (MANUAL_TRACKING==2)
                                    quedan=1;quedadist=mantrack_dist;
                                else
                                    quedan=0;quedadist=0;
                                    while(1)
                                        if (rem(point+quedan,everyNpoint)==0)
                                            quedadist=quedadist+deltaR+everyIncDec;
                                        else
                                            quedadist=quedadist+deltaR;
                                        end
                                        if (quedadist>mantrack_dist)
                                            break; % at this point, 'quedan' is the number of points
                                        else
                                            quedan=quedan+1;
                                        end
                                        
                                    end
                                end % if manual_tracking=1
                                mantrack_next_points = quedan - 1 ;
                                if (mantrack_next_points>-1)
                                    break;
                                end
                            end
                            close(h_click);
                        end % if haveImageStitched
                        
                        figure(h_focus)% gives focus to another windows because, otherwise, it seems nothing happened
                        if(finish)
                            displaysave('Breaking from mantracking-preview-window');
                            break; % OJO!! cambiado 20enero2015 para ver si así termina bien
                        end
                        if(skipping)
                            displaysave('Skipping from mantracking-preview-window');
                            break; % OJO!! cambiado 20enero2015 para ver si así termina bien
                        end
                        
                        displaysave(['manual tracking, advancing ',num2str(mantrack_dist),' mm = ', num2str(mantrack_next_points+1),' points.']);
                    end % use_recorded_track
                    
                else % the next point is already calculated
                    mantrack_next_points = mantrack_next_points - 1;
                end
                
                if (rem(point,everyNpoint)==0) % this is the nth point with possible different displacement
                    actual_deltaR = deltaR + everyIncDec;
                else
                    actual_deltaR = deltaR;
                end
                if (MANUAL_TRACKING==2)
                    actual_deltaR = mantrack_dist;
                end
                total_distance = total_distance + actual_deltaR;
                
                if (MANTRACK_RADIAL) % the X movement is actually with the angle motor
                    if (USE_RECORDED_PATH)
                        moveString = sprintf(' %0.3f %0.3f %0.3f r \r',old_recorded_path(point,1),old_recorded_path(point,2),old_recorded_path(point,3));
                    else
                        moveString = sprintf(' %0.3f %0.3f %0.3f r \r',-actual_deltaR*sin(mantrack_arctan)*dirY, 0, deltaZ);
                    end
                    fprintf(xyz,moveString); %move XYZ motor
                    if (USE_RECORDED_PATH)
                        actual_diameter = DIAMETER - param1 / param1To * (DIAMETER - DIAMETER_END );
                        ts = round( old_recorded_path(point,5) / (3.14159*actual_diameter / 2400));
                    else
                        actual_diameter = DIAMETER - param1 / param1To * (DIAMETER - DIAMETER_END ) ;
                        ts = round( -actual_deltaR*cos(mantrack_arctan)*dirX / (3.14159* actual_diameter / 2400)); % 0.15º is the single steps of the Luis' motor, 2400 steps per revolution
                    end
                    for anglestep=1:abs(ts)
                        if (ts > 0)
                            %fwrite(theta,'w');pause(0.1); total_ws_steps = total_ws_steps  + 1;
                        else
                            %fwrite(theta,'s');pause(0.1); total_ws_steps  = total_ws_steps -1;
                        end
                    end
                    if (USE_RECORDED_PATH)
                        recorded_path(point,:)=old_recorded_path(point,:);recorded_path(point,4) = defocus_mm;
                    else
                        recorded_path(point,1) = -actual_deltaR*sin(mantrack_arctan)*dirY;recorded_path(point,2) = 0; recorded_path(point,3) = deltaZ;recorded_path(point,4) = defocus_mm;recorded_path(point,5) =  -actual_deltaR*cos(mantrack_arctan)*dirX;
                    end
                    
                else % manual tracking linear 
                    if (USE_RECORDED_PATH)
                        moveString = sprintf(' %0.3f %0.3f %0.3f r \r',old_recorded_path(point,1),old_recorded_path(point,2),old_recorded_path(point,3));
                        recorded_path(point,:)=old_recorded_path(point,:);recorded_path(point,4) = defocus_mm;
                    else
                        moveString = sprintf(' %0.3f %0.3f %0.3f r \r',-actual_deltaR*sin(mantrack_arctan)*dirY, actual_deltaR*cos(mantrack_arctan)*dirX, 0);
                        recorded_path(point,1) = -actual_deltaR*sin(mantrack_arctan)*dirY;recorded_path(point,2) = actual_deltaR*cos(mantrack_arctan)*dirX;recorded_path(point,3) = deltaZ;recorded_path(point,4) = defocus_mm;recorded_path(point,5) =  0;
                    end
                    fprintf(xyz,moveString); %move XYZ motor
                    justAbsMovDone=tic(); % to check if motors are effectively stopped
                    % en modo ONLY_RECORDING_PATH pierde órdenes porque va muy rápido
                    statusXYZ=1; % v64 code to wait for motors to stop, already tested
                    while statusXYZ==1        %is moving,
                        pause(0.05);
                        if(xyz.BytesAvailable>0);fread(xyz,xyz.BytesAvailable);end
                        fprintf(xyz,' status \r');
                        statusXYZ = fscanf(xyz,'%f'); % status word
                    end
                end
            end
            if MANUAL_TRACKING
                displaysave(['  total distance= ', num2str(total_distance)]);
            end
            %%ends  code to move motors, now it is duplicated to avoid nested
        end % deltaT
        
        if(finish)
            displaysave('Breaking from offline if');
            %break OJO!! para probar
        end
        if(skipping)
            displaysave('Skipping from offline if');
            %break
        end
        
        %% ask for a direction to move over the video image
        % aqui estaba el modo manual
        
        
        point = point + 1; % incremental file naming or possible non-monotonic parameter sweep from a matrix of parameters
        
        
        if(KEYPRESS_TO_CONTINUE)
            disp('This line could set the execution into debug mode (dbstop at 918) and continue (dbcont)');
            disp('Press any key to continue, laser armed?');
            fprintf(xyz,' 1 joystick \n'); % disable joystick  EY!!! enable, just to try if it works during preview
            pause;
        end
        
        if (USE_RECORDED_PATH && point>size(old_recorded_path,1))
            % That is the last recorded point in the recorded path
            displaysave('Finish at the last recorded point');
            finish=1; % needed to avoid motor processing in param2
            break;
        end
        
        if(finish)
            break
        end
        if(skipping)
            break
        end
        
        %new v58 remaining time
        try
            now2=toc(timer_experiment); % was set at the begining
            %points_uptonow = (param1Index-1) +(param2Index-1)*Nparam1; % number of points so far
            %remainingSeconds = now2 * (Nparam1*Nparam2 - points_uptonow)/(points_uptonow);
            remainingSeconds = now2 * (effectivePoints-pointCounterEffective+1) / pointCounterEffective;
            displaysave([ '  ' datestr(datetime) ' - Remaining time (dd:hh:mm:ss): '  datestr(seconds(remainingSeconds),'DD:HH:MM:SS')]);
        catch
            displaysave(['  Remaining time: unknown']);
        end
        
        
    end % param1
    % at this point the entire param1 loop is finished, in SPLIT_PARAM2 mode we should save the spectra array before overwritten
    if (SPLIT_PARAM2==1)
        fn_split = strcat('spectra_idx2=',sprintf('%04i',param2IndexReal),'_', param2Name,'=',sprintf('%0.2f',param2));
        fn_split = strrep(fn_split,'.',','); %replaces dots with commas, just in case
        displaysave(fn_split);
        try
            save (strcat(folderPath,experimentPath,datePath,fn_split),'spectra'); % save the partial spectra for this param2 iteration
        catch
            pause(2);
            %just try again
            save (strcat(folderPath,experimentPath,datePath,fn_split),'spectra'); % save the partial spectra for this param2 iteration
        end
    end
    
    
    if(finish)
        displaysave('Breaking from param1');
        break
    end
    
    if(skipping)
        displaysave('Skipping to next param2...');
        skipping=0;
    end
    
    
    % movements AFTER each point, copied from inner loop param1
    % transverse movement is 2D scanning is activated (param1Name = deltaT)
    if ( strcmp(param2Name,'deltaT') == 1)
        % not clear what to do in this configuration
    elseif ( strcmp(param2Name,'pulsos') == 1) % en este modo se sigue capturando espectros sin moverse
        % again, if param2Name==pulsos ¿what to do?
    else
        %"normal" movement (automatic and manual) moved here
        
        %%start code to move motors, now it is duplicated to avoid nested
        %%functions
        if (MOTORS && ~MANUAL_TRACKING  && ~skipThisPoint)
            % movement
            if(USE_RECORDED_PATH)
                deltaX=-old_recorded_path(point,2);deltaY=-old_recorded_path(point,1);deltaZ=old_recorded_path(point,3);
            end
            moveString = sprintf(' %0.3f %0.3f %0.3f r \r',deltaY, deltaX, deltaZ); % EY!! OJO!! coordiantes changed from physical to image
            fprintf(xyz,moveString);
            justAbsMovDone=tic(); % to check if motors are effectively stopped
            statusXYZ=1; % v64 code to wait for motors to stop, already tested
            while statusXYZ==1        %is moving,
                pause(0.05);
                if(xyz.BytesAvailable>0);fread(xyz,xyz.BytesAvailable);end
                fprintf(xyz,' status \r');
                statusXYZ = fscanf(xyz,'%f'); % status word
            end
            actual_deltaP = deltaP;
            if (deltaP ~= 0.0) % angle motor
                if (rem(point,everyNpoint)==0) % this is the nth point with possible different displacement
                    actual_deltaP = deltaP + everyIncDec;
                end
                if (USE_RECORDED_PATH)
                    actual_deltaP = old_recorded_path(point,5);
                end
                total_distance = total_distance + actual_deltaP + sqrt(deltaX^2+deltaY^2);
                actual_diameter = DIAMETER - param1 / param1To * (DIAMETER - DIAMETER_END );
                theta_steps = round( actual_deltaP / (3.14159*actual_diameter / 2400));
                for anglestep=1:theta_steps
                    if (deltaP > 0)
                        %fwrite(theta,'w');pause(0.1); total_ws_steps = total_ws_steps + 1;
                    else
                        %fwrite(theta,'s');pause(0.1); total_ws_steps = total_ws_steps - 1 ;
                    end
                end
            end
            recorded_path(point,1) = deltaX;recorded_path(point,2) = deltaY;recorded_path(point,3) = deltaZ;recorded_path(point,4) = defocus_mm;recorded_path(point,5) = actual_deltaP;
        end
        
        if (MOTORS && MANUAL_TRACKING)
            if(mantrack_next_points == 0)  % there are no automatic points to follow the path
                if (~USE_RECORDED_PATH) % not using using_recorded_path
                    sound(sound2);sound(sound3); % sound alert
                    if (0)  % PREVIEW shows the video preview to adjust the focus is needed
                        h_preview=figure('Position',[100 100 640  480],'Name','Video Preview'); % create preview window
                        hold on;
                        displaysave('Showing video preview (flipped), press any to close and choose NEW DIRECTION...');
                        k=[]; set(h_preview,'keypress','k=get(gcf,''currentchar'');');
                        previewTime=tic;
                        frame = captureFrameGuppy;
                        h_ims_prev=imshow(frame,'InitialMagnification', 67);
                        while 1
                            frame = captureFrameGuppy;
                            set(h_ims_prev,'CData',frame); % new! v28 update image without a new imshow, promise speed things up AND avoid memory leak error
                            if (0 && toc(previewTime) > 300) % leak memory problem
                                sound(sound2);sound(sound3); % sound alert
                                disp('Automatic pause to prevent memory leak error, press any key to continue...');
                                pause;
                                previewTime=tic;
                            end
                            %set(h_preview,'ButtonDownFcn','disp(''ay!'')');
                            target(pixel_shotX,pixel_shotY); % draw a cross around
                            if (ESCALA==1) % draw 1mm ticks starting at the laser spot
                                if calibrated_pixels_mm_XY_X > 300   % new v58
                                    mm1o01=calibrated_pixels_mm_XY_X*0.1; % line every 100um
                                    numberOflines=7;
                                else
                                    mm1o01=calibrated_pixels_mm_XY_X; % line every 1mm
                                    numberOflines=3;
                                end
                                hold on;
                                for horticks=1:numberOflines % cortos hacia la izquierda
                                    line([pixel_shotX-(horticks-1)*mm1o01  pixel_shotX-(horticks-1)*mm1o01],[0 300],'Color','w');
                                    line([pixel_shotX-(horticks-1)*mm1o01+1  pixel_shotX-(horticks-1)*mm1o01+1],[0 300],'Color','k');
                                end
                                for horticks=1:numberOflines  % largos hacia la derecha
                                    line([pixel_shotX+(horticks)*mm1o01  pixel_shotX+(horticks)*mm1o01],[0 600],'Color','w');
                                    line([pixel_shotX+(horticks)*mm1o01+1  pixel_shotX+(horticks)*mm1o01+1],[0 600],'Color','k');
                                end
                            end
                            if (~isempty(k) || ~ishandle(h_preview))  % if it has been closed
                                break;
                            end
                        end
                        close(h_preview);
                    end % preview
                    
                    %twitter
                    if (TWITTER==2)
                        try
                            tw.updateStatus(['Autolibs: requiere nueva dirección manual_tracking ' datestr(now) ]);
                        catch
                            TWITTER=0;
                        end
                    end
                    
                    h_click=figure('Position',[1 1 scrsz(3) scrsz(4)],'Name','Current image - click for direction of new points'); % create new large window
                    pause(1); % new method (v27) return the LAST capture frame, it could be captured WHILE the motors are moving
                    frame = captureFrameGuppy;
                    imshow(frame,'InitialMagnification', 67); hold on;
                    target(pixel_shotX,pixel_shotY); % draw a cross around
                    if (ESCALA==1) % draw 1mm ticks starting at the laser spot
                        if calibrated_pixels_mm_XY_X > 300   % new v58
                            mm1o01=calibrated_pixels_mm_XY_X*0.1; % line every 100um
                            numberOflines=7;
                        else
                            mm1o01=calibrated_pixels_mm_XY_X; % line every 1mm
                            numberOflines=3;
                        end
                        hold on;
                        for horticks=1:numberOflines % cortos hacia la izquierda
                            line([pixel_shotX-(horticks-1)*mm1o01  pixel_shotX-(horticks-1)*mm1o01],[0 300],'Color','w');
                            line([pixel_shotX-(horticks-1)*mm1o01+1  pixel_shotX-(horticks-1)*mm1o01+1],[0 300],'Color','k');
                        end
                        for horticks=1:numberOflines  % largos hacia la derecha
                            line([pixel_shotX+(horticks)*mm1o01  pixel_shotX+(horticks)*mm1o01],[0 600],'Color','w');
                            line([pixel_shotX+(horticks)*mm1o01+1  pixel_shotX+(horticks)*mm1o01+1],[0 600],'Color','k');
                        end
                    end
                    
                    while(1)
                        [xc, yc] = ginput(1);
                        mantrack_dist =  sqrt(((xc-pixel_shotX)/calibrated_pixels_mm_XY_X)^2+((yc-pixel_shotY)/calibrated_pixels_mm_XY_Y)^2); % distance in mm
                        mantrack_arctan = atan(abs(yc-pixel_shotY)/abs(xc-pixel_shotX)); % tangent(angle)
                        if (xc>pixel_shotX)
                            dirX=1;
                        else
                            dirX=-1;
                        end
                        if (yc>pixel_shotY)
                            dirY=1;
                        else
                            dirY=-1;
                        end
                        % wee need to calculate how many points there is until next
                        % stop, taking into acount everyNpoints;
                        quedan=0;quedadist=0;
                        while(1)
                            if (rem(point+quedan,everyNpoint)==0)
                                quedadist=quedadist+deltaR+everyIncDec;
                            else
                                quedadist=quedadist+deltaR;
                            end
                            if (quedadist>mantrack_dist)
                                break; % at this point, 'quedan' is the number of points
                            else
                                quedan=quedan+1;
                            end
                            
                        end
                        mantrack_next_points = quedan - 1 ;
                        if (mantrack_next_points>-1)
                            break;
                        end
                    end
                    close(h_click);
                    
                    if(finish)
                        displaysave('Breaking from use-recorded-track if');
                        break;
                    end
                    
                    displaysave(['manual tracking, advancing ',num2str(mantrack_dist),' mm = ', num2str(mantrack_next_points+1),' points.']);
                end % use_recorded_track
                
            else % the next point si already calculated
                mantrack_next_points = mantrack_next_points - 1;
            end
            
            if (rem(point,everyNpoint)==0) % this is the nth point with possible different displacement
                actual_deltaR = deltaR + everyIncDec;
            else
                actual_deltaR = deltaR;
            end
            
            total_distance = total_distance + actual_deltaR;
            
            if (MANTRACK_RADIAL) % the X movement is actually with the angle motor
                if (USE_RECORDED_PATH)
                    moveString = sprintf(' %0.3f %0.3f %0.3f r \r',old_recorded_path(point,1),old_recorded_path(point,2),old_recorded_path(point,3));
                else
                    moveString = sprintf(' %0.3f %0.3f %0.3f r \r',-actual_deltaR*sin(mantrack_arctan)*dirY, 0, deltaZ);
                end
                fprintf(xyz,moveString); %move XYZ motor
                if (USE_RECORDED_PATH)
                    actual_diameter = DIAMETER - param1 / param1To * (DIAMETER - DIAMETER_END );
                    ts = round( old_recorded_path(point,5) / (3.14159*actual_diameter / 2400));
                else
                    actual_diameter = DIAMETER - param1 / param1To * (DIAMETER - DIAMETER_END ) ;
                    ts = round( -actual_deltaR*cos(mantrack_arctan)*dirX / (3.14159* actual_diameter / 2400)); % 0.15º is the single steps of the Luis' motor, 2400 steps per revolution
                end
                for anglestep=1:abs(ts)
                    if (ts > 0)
                        %fwrite(theta,'w');pause(0.1); total_ws_steps = total_ws_steps  + 1;
                    else
                        %fwrite(theta,'s');pause(0.1); total_ws_steps  = total_ws_steps -1;
                    end
                end
                if (USE_RECORDED_PATH)
                    recorded_path(point,:)=old_recorded_path(point,:);recorded_path(point,4) = defocus_mm;
                else
                    recorded_path(point,1) = -actual_deltaR*sin(mantrack_arctan)*dirY;recorded_path(point,2) = 0; recorded_path(point,3) = deltaZ;recorded_path(point,4) = defocus_mm;recorded_path(point,5) =  -actual_deltaR*cos(mantrack_arctan)*dirX;
                end
                
            else % manual tracking linear
                if (USE_RECORDED_PATH)
                    moveString = sprintf(' %0.3f %0.3f %0.3f r \r',old_recorded_path(point,1),old_recorded_path(point,2),old_recorded_path(point,3));
                    recorded_path(point,:)=old_recorded_path(point,:);recorded_path(point,4) = defocus_mm;
                else
                    moveString = sprintf(' %0.3f %0.3f %0.3f r \r',-actual_deltaR*sin(mantrack_arctan)*dirY, actual_deltaR*cos(mantrack_arctan)*dirX, 0);
                    recorded_path(point,1) = -actual_deltaR*sin(mantrack_arctan)*dirY;recorded_path(point,2) = actual_deltaR*cos(mantrack_arctan)*dirX;recorded_path(point,3) = deltaZ;recorded_path(point,4) = defocus_mm;recorded_path(point,5) =  0;
                end
                fprintf(xyz,moveString); %move XYZ motor
            end
        end
        
    end % deltaT
    
    
end % param2


%% experiment done
%COPIAR DESDE AQUÍ HASTA EL FINAL SI EL PROGRAMA NO TERMINA
% save the spectral window information and other
ancho = lambdas(Npixeles)-lambdas(1);
central = ceil(lambdas(1) + ancho/2);
nm2pixels=1024/ancho;

if (ancho > 250) % grating 150lpp
    if (USE_AVANTES)
        grating=['avantes 8ch'];
    else
        grating=['150lpp@' num2str(central) 'nm'];
    end
else
    if (ancho > 30) % graging 1200lpp
        grating=['1200lpp@' num2str(central) 'nm'];
    else
        grating=['1800lpp@' num2str(central) 'nm'];
    end
end

%displaysave(['\nWindow ' grating '\n']);
%displaysave(['first lambda ' num2str(lambdas(1)) '\n']);
%displaysave(['last  lambda ' num2str(lambdas(Npixeles))  '\n']);
%displaysave(['Gain  x' num2str(gainIntensifier)  '\n']);

if (finish && ~panicAvantes)
    displaysave('Finished by user');
end

if (finish && panicAvantes)
    displaysave('Finish due to Avantes panic');
end

if(USE_OCEANOPTICS)
    wrapper.closeAllSpectrometers();
end

% save matlab for data for further processing
saveas(h_focus,strcat(folderPath,experimentPath,datePath,'ventana_ratio.png'),'png');
close(h_focus);  % or don't close
stop(vid);
if (TAKE_WEBCAM_PHOTO==1)
    stop(vid2);
    clear vid2;
end
clear vid; % this variable gives warnings when loading the file
clear src;
clear camaras;
clear h_danger;

displaysave('Saving everything, please wait...');

try % window could be already closed
    if havePolygon % save roi variable just in case
        save (strcat(folderPath,experimentPath,datePath,'polygon.mat'),'roiPol');
    end
    if haveLinePath
        save (strcat(folderPath,experimentPath,datePath,'linePath.mat'),'roiLine');
    end
    if haveImageStitched
        savefig(h_stitch,strcat(folderPath,experimentPath,datePath,'stitchedImage.fig'));
        imwrite(stitchImage, strcat(folderPath,experimentPath,datePath,'stitchedImage.png')); % full resolution image
    end
catch
end

clear h_stitch; % new v60 got saved with warning due to size

%clear wrapper; % no se puede borrar porque ya no vuelve a funcionar!!!
%master save of matlabdata new v63 do not save stitchimage, it could be very large is already saved in a .fig file
save (strcat(folderPath,experimentPath,datePath,dataFile),'-regexp','^(?!(wrapper|vid*|obj*|stitchImage)$).','-v7.3');
save (strcat(folderPath,experimentPath,datePath,'recorded_path'),'recorded_path'); % save the recorded path so it can be used again if needed


%% CIERRE FINAL
if (USE_ARDUINO_QSWITCH)
    fprintf(qswitch,'0');  % cut the q-switch signal
    fclose(qswitch);
    if (USE_DUAL_PULSE_LASER == 0) % for DP laser shutter and qswitch is controlled by software
        h=warndlg('¡ATENCIÓN! apagar el shutter y el Q-switch en el mando del láser.');
    end
    %uiwait(h);
end
if (USE_PULSE_GENERATOR==1) fprintf(pulsegen,['q0' char(13)]); fclose(pulsegen); end % activates qswitch pulses with no 'noise' pulse

if (USE_DUAL_PULSE_LASER)
    if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
    pause(0.1);
    fwrite(lotis,['STOP' char(13)]);    % new v55 stop the lamp, so the experiments could be left running
    pause(0.1);
    if(lotis.BytesAvailable>0);fread(lotis,lotis.BytesAvailable);end
    fclose(lotis);
end
if MOTORS
    moveString = sprintf('%f %f %f move \n',backToStartingCoordinates(1),backToStartingCoordinates(2),backToStartingCoordinates(3)); % new v62 back to where it was, home is where it was
    fprintf(xyz,moveString); pause(0.05);
    fprintf(xyz,' 1 joystick \n'); % enable joystick
    total_ws_steps = rem(total_ws_steps,2400); % remove entire turns
    if (total_ws_steps > 30)
        disp('Be patient, there are many motor steps to bring it back to the origin...');
    end
    for anglestep=1:abs(total_ws_steps) % send stepper motor back to the original position
        if (total_ws_steps < 0)
            %fwrite(theta,'w'); pause(0.03);
        else
            %fwrite(theta,'s'); pause(0.03);
        end
    end
    pause(0.1);
    fprintf(xyz,' 1 joystick \n'); % enable joystick, again, just in case, it fails from time to time
    %fclose(theta);
    fclose(xyz);
end

% new v64 to check movements are ok
if (~ ONLY_RECORDING_PATH && CAMERA) % if only_recording_path no photos to speed up things
    pause(3); % going back can take a while
    frameAfter = captureFrameGuppy;
    imwrite(frameAfter,strcat(folderPath,experimentPath,datePath,'AfterPhoto.png')); %save current frame with the same name as the experiment
end

if (USE_LEDLIGHT)
    fclose(ledlight);
end
if (USE_SHUTTER_CTRL)
    fprintf(shutter,'ens?');pause(0.1);if(shutter.BytesAvailable>0);response=fread(shutter,shutter.BytesAvailable);end
    if (numel(response) >= 6 && response(6)=='1')
        fprintf(shutter,'ens');pause(0.1);if(shutter.BytesAvailable>0);response=fread(shutter,shutter.BytesAvailable);end % open the shutter
    end
    fclose(shutter);
end

if (USE_AVANTES)
    spectrometer('deactivate',h1);
    spectrometer('deactivate',h2);
    spectrometer('deactivate',h3);
    spectrometer('deactivate',h4);
    spectrometer('deactivate',h5);
    spectrometer('deactivate',h6);
    spectrometer('deactivate',h7);
    spectrometer('deactivate',h8);
    spectrometer('done');
end
% close danger window, we always forget to do it
try
    close(h_danger)
catch
end

%twitter
if (TWITTER>0)
    try
        tw.updateStatus(['Autolibs: terminado ' datestr(now) ]);
    catch
    end
end

if (~ ONLY_RECORDING_PATH) % if only_recording_path, there is no lambdas vectors and procesaMAT crashes
    DESDE_AUTOLIBS = 1;
end

displaysave(['Estimated time was: ' datestr(seconds(timePerPoint*effectivePoints),'DD:HH:MM:SS')]); % if effectivePoints is not defined, this is because is only calculated for PosX/PosY, it should be calc from number of points from the idx1/idx2 loops, not done yet
totalSeconds = 86400*(now-whenStarted);
displaysave(['Total time was: ' datestr(totalSeconds/86400,'DD:HH:MM:SS')]);
displaysave(['Overhead time was: ' num2str(totalSeconds/effectivePoints - 0.1 * Npulses) ]);

try
    tgprintf([ 'Terminada medida LIBS ...' SUBEXPERIMENT_NAME ' en ' getenv('COMPUTERNAME')]);
catch
end

disp('Done!');

% new v65 function that skip points inside and outside a 3mm ring
function skip = ring3mm(x,y,p1f,p1t,p1s,p2f,p2t,p2s)
% loop should be a 3mm square
x0 = (p1f+p1t)/2; % center of square
y0 = (p2f+p2t)/2; 
r =sqrt((x-x0)^2 + (y-y0)^2) ; 
if r < 1 || r > 1.5
    skip = 1;
else
    skip = 0;
end
end