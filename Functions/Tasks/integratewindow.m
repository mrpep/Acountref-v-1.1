function [integratedir] = integratewindow(ir,fs,tau)

    begink = FindBeginIR(ir,0);
    ir = ir(begink:end);
    
%Ventana:
    recttime = 0.01; %Tiempo de ventana rectangular
    %Sección rectangular:
    rectwin = ones(1,fs*recttime);
    %Sección exponencial:
    expwin = zeros(1,floor(5*tau*fs));
    for i= 1 : floor(5*tau*fs);
       expwin(i) = exp(-i/floor(tau*fs));
    end 
    window = horzcat(rectwin,expwin);
    %windowzeros = zeros(1,2*length(window));
    %length(ir)
    %ir = horzcat(ir,windowzeros);
    %length(ir)
    %Division en parte positiva y negativa:
    ir_positive=ir;
    ir_negative=ir;
    ir_positive(ir<0)=0;
    ir_negative(ir>0)=0;
    
    %ir_negative=(ir_negative).*(-1);
    %window = window';
    %Integración:
    for i = 1:length(ir)-length(window)
        ir_positive_windowed=(ir_positive(i:i+length(window)-1)).*window;
        ir_negative_windowed=-(((ir_negative(i:i+length(window)-1)).*window));           
            
        
        ir_positive_windowedS(i)=sum(ir_positive_windowed);
        ir_negative_windowedS(i)=sum(ir_negative_windowed);
        
    end       

    
    ir_negative_windowedS = -ir_negative_windowedS;
    integratedir=(ir_positive_windowedS)'+(ir_negative_windowedS)';
    %zerostoadd = zeros(irlength - length(integratedir),1);

    %integratedir = vertcat(integratedir,zerostoadd);
    integratedir = integratedir';   
    %integratedir = integratedir(1:end-length(window));
    %integratedir = integratedir - mean(integratedir);
    %length(ir)
    integratedir = integratedir - mean(integratedir);
    integratedir=integratedir./max(abs(integratedir));
    

end