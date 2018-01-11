clear all
clc



%% Function call
clearvars;close all;clc;
rng(1) ;% to ensure reproducibility of the example.
filename        = 'INPUT.txt';
[u,v,w,t,nodes] = windSim(filename);
U               = u+nodes.U*ones(size(t)); % get fluctuating + mean value
vel             = u; %REVISAR EN LA FUNCIÓN WINDSIM QUE ESPECTRO DE DENSIDAD DE POTENCIA DEL VIENTO SE UTILIZÓ PARA QUE COINCIDA CON LA OPCIÓN ELEGIDA MANUALMENTE EN ESTE PROGRAMA



%% Datos
At     = [640,640,640,640,640,640,640,640,640,640];    %Áreas tributarias para cada entrepiso

zo     = 1;                   
zr     = 10;
T      = 600;                              
U_m_zr = 30;                       
pa     = 1.225;                           
CD     = 0.8;

% Resultados y comparaciones
z    = nodes.Z;
U_m  = nodes.U;
alfa = 0.096*(log10(zo))+0.016*((log10(zo))^2)+0.24;
for i = 1:length(z)
    if z(i)<10
       U_m(i)=U_m_zr;
    else
       U_m(i)=U_m_zr*((z(i)/zr)^alfa);
    end
end
v = 0.67+0.05*log(zo);
for i = 1:length(z)
   u_f(i) = (0.4*U_m(i))/(log(z(i)/zo));   %Velocidad de fricción
   de(i)  = (5.7*(u_f(i)^2))^(1/2);        %Desviación estandar
   Iu(i)  = de(i)/U_m(i);                  %Intensidad de turbulencia
   Lu(i)  = 300*((z(i)/200)^v);            %Longitud de turbulencia
end

f1     = 0.001;
fu     = 10;
deltat = 1/(2*fu);
fs     = 1/deltat;
N      = (T/deltat)+1;                         % Número de datos de las simulaciones
deltaf = (fu-f1)/(N-1);
deltaw = 2*pi*deltaf;
f(1)   = f1;
for i = 2:N
    f(i) = f(i-1)+deltaf;
end
w = f.*2*pi;
t(1) = 0;
t(2) = deltat;
for i = 3:N
    t(i) = t(i-1)+deltat;
end

fprintf('\n');
disp('Elija la función de densidad espectral de potencia a utilizar:');
fprintf('\n');
disp('Opción 1: Von Kármán');
disp('Opción 2: Von Kárman-Harris');
disp('Opción 3: Kaimal');
disp('Opción 4: Kaimal modificado');
disp('Opción 5: Davenport');
disp('Opción 6: Solari');
fprintf('\n');
opcion = input('Presione 1,2,3,4,5 ó 6 según sea la opción a elegir:');
fprintf('\n');

%Función de densidad de potencia espectral del viento
%Se hace a lo alto del edificio (z) y a las diferentes frecuencias (N)
Su = zeros(length(z),N);
if opcion == 1
    for i = 1:length(z)
        for j = 1:N
            Su(i,j) = ((4*5.7*(u_f(i)^2)*Lu(i))/U_m(i))/(1.339*((1+(39.48*(((f(j)*Lu(i))/U_m(i))^2)))^(5/6)));
        end
    end
    
elseif opcion == 2
    for i = 1:length(z)
        for j = 1:N
            Su(i,j) = ((4*(de(i)^2)*Lu(i))/U_m(i))/(1*((1+(70.8*(((f(j)*Lu(i))/U_m(i))^2)))^(5/6)));
        end
    end
    
elseif opcion == 3
    for i = 1:length(z)
        for j = 1:N
            Su(i,j) = ((105*(u_f(i)^2)*z(i))/U_m(i))/(1*((1+(33*(((f(j)*z(i))/U_m(i))^1)))^(5/3)));
        end
    end
  
elseif opcion == 4
    for i = 1:length(z)
        for j = 1:N
            Su(i,j) = ((100*(u_f(i)^2)*z(i))/U_m(i))/(1*((0.44+(33*(((f(j)*z(i))/U_m(i))^1)))^(5/3)));
        end
    end
    
elseif opcion == 5
    for i = 1:length(z)
        for j = 1:N
            Su(i,j) = (2*(((f(j)*1200)/U_m(i))^2)*(de(i)^2))/(f(j)*(3*((1+(((f(j)*1200)/U_m(i))^2))^(4/3))));
        end
    end
    
elseif opcion == 6
    for i = 1:length(z)
        for j = 1:N
            Su(i,j) = ((6.868*(de(i)^2)*Lu(i))/U_m(i))/(1*((1+(10.302*(((f(j)*Lu(i))/U_m(i))^1)))^(5/3)));
        end
    end
    
end

fprintf('\n')
disp('Las desviaciones estándar buscadas son:');
de = de'
fprintf('\n')
disp('Las velocidades medias del viento buscadas son:');
U_m
fprintf('\n')
disp('Las intensidades de turbulencia buscadas son:');
Iu = Iu'
fprintf('\n')
disp('Las longitudes de turbulencia buscadas son:');
Lu = Lu'

%Desviación estándar del viento a lo alto del edificio (z) y a diferentes frecuencias (N) 
fprintf('\n')
disp('Las desviaciones estándar obtenidas de las simulaciones son:');
for j = 1:length(z)
    for i = 1:N
        suma(i) = (vel(j,i)^2);
    end
    de_o(j) = ((1/(N-1))*sum(suma))^(1/2);
end
de_o = de_o'


%Función de admitancia aerodinámica a lo alto a diferentes frecuencias
adm = zeros(length(z),N);                           %Función de admitancia aerodinámica
for j = 1:length(z)
    for i = 1:N
        adm(j,i) = (1+(((2*f(i)*(At(j)^(1/2)))/U_m(j))^(4/3)))^(-7/6);
    end
end

Su_real   = zeros(length(z),N);         % Transformada rápida de Fourier obtenida del registro de velocidad 
Su_real_f = zeros(length(z),N);         % Transformada rápida de Fourier obtenida del registro de velocidad filtrado 
vel_f     = zeros(length(z),N);         %Registro de velocidad del viento filtrado mediante la función de admitancia aerodinámica 
for i=1:length(z)
    Su_real(i,:)   = fft(vel(i,:)); 
    Su_real_f(i,:) = Su_real(i,:).*adm(i,:);
    vel_f(i,:)     = ifft(Su_real_f(i,:)); 
end


vel_t = zeros(length(z),N);                 %Registro de velocidades totales del viento sin filtro
for i = 1:length(z)
    vel_t(i,:) = vel(i,:)+U_m(i); 
end

disp('Las velocidades media obtenidas de las simulaciones son:');
for i = 1:length(z)
    U_m_o(i) = (1/N)*sum(vel_t(i,:));
end
U_m_o = U_m_o'

fprintf('\n')
disp('Las intensidades de turbulencia obtenidas de las simulaciones son:');
for i = 1:length(z)
    Iu_o(i) = de_o(i)/U_m_o(i);
end
Iu_o = Iu_o'

vel_t_f = zeros(length(z),N);               %Registro de velocidades totales del viento con filtro
for i = 1:length(z)
    vel_t_f(i,:) = vel_f(i,:)+U_m(i); 
end

corr = zeros(length(z),N);                 %Función de autocorrelación completa
for i = 1:length(z)
    corr(i,:) = autocorr(vel_t(i,:),(N-1));
end

corre   = zeros(length(z),N);                %Función de autocorrelación hasta cortar en cero
retraso = zeros(length(z),N);
for j = 1:length(z)
    for i = 1:N
        if corr(j,i) > 0
           corre(j,i)   = corr(j,i);                
           retraso(j,i) = t(i);
           continue
        else
           break
        end
    end
end

for i = 1:length(z)                       %Escala de tiempo
    Tu(i) = trapz(retraso(i,:),corre(i,:));
end

fprintf('\n')
disp('Las longitudes de turbulencia obtenidas de la simulaciones son:');
for i = 1:length(z)
    Lu_o(i) = Tu(i)*U_m_o(i);               
end
Lu_o = Lu_o'

% Espectros de densidad de potencia de las simulaciones

for i = 1:length(z)
   [PSD(i,:),ff(i,:)] = pwelch(vel(i,:),N,round(N/2),N,fs);       %Sin filtro
end

for i = 1:length(z)
   [PSD_f(i,:),fff(i,:)] = pwelch(vel_f(i,:),N,round(N/2),N,fs);   %Con filtro
end

%Gráficas

fprintf('\n');
disp('¿ Desea graficar algún registro de velocidades del viento ?:');
fprintf('\n');
disp('Opción 1: Si');
disp('Opción 2: No');
fprintf('\n');
opc = input('Presione 1 ó 2 según sea la opción a elegir:');
fprintf('\n');
if opc == 1
    
    altura = input('Ingrese la altura medida desde el suelo donde quiere graficar:');
    for i = 1:length(z)
        if altura == z(i)
            altura = i;
            break
        else
            continue
        end
    end
    
    subplot(3,1,2)    
    plot(f,adm(altura,:),'g-')  
    xlabel('\it f (Hz)'); 
    ylabel('\it X (z,f)^2');
    title('Función de admitancia aerodinámica');
    grid on
    set(gcf,'Color','w')
    
    subplot(3,1,1)    
    plot(t,vel_t(altura,:),'r-')  
    xlabel('\it t (s)'); 
    ylabel('\it U (m/s)');
    title('Velocidades longitudinales del viento');
    grid on
    set(gcf,'Color','w')
    
    subplot(3,1,3)    
    plot(t,vel_t_f(altura,:),'b-')  
    xlabel('\it t (s)'); 
    ylabel('\it U (m/s)');
    title('Velocidades longitudinales filtradas del viento');
    grid on
    set(gcf,'Color','w') 
else
end

fprintf('\n');
disp('¿ Desea graficar el espectro de amplitudes de algún registro de velocidades ?:');
fprintf('\n');
disp('Opción 1: Si');
disp('Opción 2: No');
fprintf('\n');
op=input('Presione 1 ó 2 según sea la opción a elegir:');
fprintf('\n');
if op == 1
    
     altura = input('Ingrese la altura medida desde el suelo donde quiere graficar:');
    for i = 1:length(z)
        if altura == z(i)
            altura = i;
            break
        else
            continue
        end
    end
    
    loglog(ff(altura,:),PSD(altura,:),'r-')  
    hold on   
    loglog(fff(altura,:),PSD_f(altura,:),'b-')
    loglog(f,Su(altura,:),'g-') 
    xlabel('\it f (Hz)'); 
    ylabel('\it Su (m^2/s^2Hz)'); 
    axis ([f1 fu  (1/10)*min(PSD_f(altura,:)) 10*max(PSD(altura,:))])
    legend('\it Sin filtro','\it Con filtro','\it Teórico');
    title('Espectro de densidad de potencia del viento');
    grid on
    set(gcf,'Color','w') 
else
end

F = zeros(length(z),N);                    %Fuerza longitudinal generada por el viento sin filtro
for i = 1:length(z)
    F(i,:) = (1/2)*pa*At(i)*CD.*(vel_t(i,:).^2);
end

F_f = zeros(length(z),N);                  %Fuerza longitudinal generada por el viento con filtro
for i = 1:length(z)
    F_f(i,:) = (1/2)*pa*At(i)*CD.*(vel_t_f(i,:).^2);
end
       

fprintf('\n');
disp('¿ Desea graficar el registro de fuerzas del viento ?:');
fprintf('\n');
disp('Opción 1: Si');
disp('Opción 2: No');
fprintf('\n');
opci = input('Presione 1 ó 2 según sea la opción a elegir:');
fprintf('\n');
if opci == 1
    
    altura = input('Ingrese la altura medida desde el suelo donde quiere graficar:');
    for i = 1:length(z)
        if altura == z(i)
            altura = i;
            break
        else
            continue
        end
    end
    
    subplot(3,1,2)    
    plot(f,adm(altura,:),'g-')  
    xlabel('\it f (Hz)'); 
    ylabel('\it X (z,f)^2');
    title('Función de admitancia aerodinámica');
    grid on
    set(gcf,'Color','w')
    
    subplot(3,1,1)    
    plot(t,F(altura,:),'r-')  
    xlabel('\it t (s)'); 
    ylabel('\it F (N)'); 
    title('Fuerzas longitudinales del viento');
    grid on
    set(gcf,'Color','w')
    
    subplot(3,1,3)    
    plot(t,F_f(altura,:),'b-')  
    xlabel('\it t (s)'); 
    ylabel('\it F (N)'); 
    title('Fuerzas longitudinales filtradas del viento');
    grid on
    set(gcf,'Color','w') 
else
end

cov    = zeros(length(z)-1,1);     % Covarianza entre entrepisos consecutivos
c_corr = zeros(length(z)-1,1);  % Coeficiente de correlación entre entrepisos consecutivos
vv     = zeros(length(z),N);
dede   = zeros(length(z)-1,1);
for i = 1:(length(z)-1)
    vv(i,:)   = vel(i,:).*vel(i+1,:);
    cov(i)    = (1/T)*trapz(t,vv(i,:));
    dede(i)   = de_o(i)*de_o(i+1);
    c_corr(i) = cov(i)/dede(i);
end
fprintf('\n')
disp('Los coeficientes de correlación entre entrepisos consecutivos son:');
c_corr

%Exportar resultados a Excel 

Matriz      = zeros(length(z),5);
Matriz(:,1) = z;
Matriz(:,2) = de;
Matriz(:,3) = U_m;
Matriz(:,4) = Iu;
Matriz(:,5) = Lu;

Matriz2      = zeros(length(z),5);
Matriz2(:,1) = z;
Matriz2(:,2) = de_o;
Matriz2(:,3) = U_m_o;
Matriz2(:,4) = Iu_o;
Matriz2(:,5) = Lu_o;

Matriz3                    = zeros(N,(length(z)+1));
Matriz3(:,1)               = f;
Matriz3(:,2:(length(z)+1)) = adm';

Matriz4                    = zeros(N,(length(z)+1));
Matriz4(:,1)               = f;
Matriz4(:,2:(length(z)+1)) = Su';

Matriz5                    = zeros(length(PSD(1,:)),(length(z)+1));
Matriz5(:,1)               = ff(1,:);
Matriz5(:,2:(length(z)+1)) = PSD';

Matriz6                    = zeros(length(PSD_f(1,:)),(length(z)+1));
Matriz6(:,1)               = fff(1,:);
Matriz6(:,2:(length(z)+1)) = PSD_f';

Matriz7                    = zeros(N,2*length(z));

j       = 1;
retraso = retraso';
for i=1:length(z)
    Matriz7(:,j) = retraso(:,i);
    j            = j+2;
end
j     = 2;
corre = corre';
for i = 1:length(z)
    Matriz7(:,j) = corre(:,i);
    j            = j+2;
end

Matriz8                    = zeros(N,(length(z)+1));
Matriz8(:,1)               = t;
Matriz8(:,2:(length(z)+1)) = vel_t';

Matriz9                    = zeros(N,(length(z)+1));
Matriz9(:,1)               = t;
Matriz9(:,2:(length(z)+1)) = vel_t_f';

Matriz10                   = zeros(N,(length(z)+1));
Matriz10(:,1)              = t;
Matriz10(:,2:(length(z)+1))= F';

Matriz11                   = zeros(N,(length(z)+1));
Matriz11(:,1)              = t;
Matriz11(:,2:(length(z)+1))=F_f';


texto1 = {'Altura','Des.est.','Vel.media','Iu','Lu'};
texto2 = {'f','X(z,f)^2'};
texto3 = {'f','Su'};
texto4 = {'Retraso','Autocorrelación','Retraso','Autocorrelación','...','...'};
texto5 = {'Coeficientes de correlación entre alturas consecutivas'};
texto6 = {'Tiempo', 'Velocidades'};
texto7 = {'Tiempo','Fuerzas'};

xlswrite('Características del viento.xlsx',texto1,'Valores esperados','A1');
xlswrite('Características del viento.xlsx',Matriz,'Valores esperados','A2');
xlswrite('Características del viento.xlsx',texto1,'Valores obtenidos','A1');
xlswrite('Características del viento.xlsx',Matriz2,'Valores obtenidos','A2');

xlswrite('Admitancia aerodinámica.xlsx',texto2,'Filtro','A1');
xlswrite('Admitancia aerodinámica.xlsx',Matriz3,'Filtro','A2');

xlswrite('Espectros de densidad de potencia del viento.xlsx',texto3,'Teórico','A1');
xlswrite('Espectros de densidad de potencia del viento.xlsx',Matriz4,'Teórico','A2');
xlswrite('Espectros de densidad de potencia del viento.xlsx',texto3,'Sin filtro','A1');
xlswrite('Espectros de densidad de potencia del viento.xlsx',Matriz5,'Sin filtro','A2');
xlswrite('Espectros de densidad de potencia del viento.xlsx',texto3,'Con filtro','A1');
xlswrite('Espectros de densidad de potencia del viento.xlsx',Matriz6,'Con filtro','A2');

xlswrite('Autocorrelación.xlsx',texto4,'Funciones','A1');
xlswrite('Autocorrelación.xlsx',Matriz7,'Funciones','A2');
xlswrite('Autocorrelación.xlsx',texto5,'Coeficientes','A1');
xlswrite('Autocorrelación.xlsx',c_corr,'Coeficientes','A2');

xlswrite('Velocidades del viento.xlsx',texto6,'Sin filtro','A1');
xlswrite('Velocidades del viento.xlsx',Matriz8,'Sin filtro','A2');
xlswrite('Velocidades del viento.xlsx',texto6,'Con filtro','A1');
xlswrite('Velocidades del viento.xlsx',Matriz9,'Con filtro','A2');

xlswrite('Fuerzas del viento.xlsx',texto7,'Sin filtro','A1');
xlswrite('Fuerzas del viento.xlsx',Matriz10,'Sin filtro','A2');
xlswrite('Fuerzas del viento.xlsx',texto7,'Con filtro','A1');
xlswrite('Fuerzas del viento.xlsx',Matriz11,'Con filtro','A2');