function [u,v,w,t,nodes] = windSim(filename)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% [u,v,w,t,nodes] = windSim(filename) generates spatially correlated wind 
	% histories based on an input file "filename" which is a text file. 
	% The space is defined using a cartesian coordinate system (x,y,z).
	% The x axis is the horizontal axis aligned with the wind direction.
	% The y axis is the horizontal axis normal to the wind direction. 
	% The z axis is the veertical axis. 
	% In this simulation the flow is assumed to have only a mean value in the x
	% direction. In other words, mean(v)=mean(w)=0 m/s. This is a common
	% assumption in wind engineering
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% INPUT: text file .
	% example: filename = 'INPUT.txt'
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% OUTPUT:
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% u: [Nyy x N] matrix of the along wind component
	% v: [Nyy x N] matrix of the cross-wind component
	% w: [Nyy x N] matrix of the vertical wind component
	% t: time vector
	% nodes: structure variables that contains informations about mean wind
	% speed, and coordinates of each nodes
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% example: [u,w,t,nodes] = windSim('INPUT.txt');
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%  Author. Etienne Cheynet  -- last modified: 22/04/2016
	%  see also windSim2.m coherence.m
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%% TIME DEFINITION
	% import data from input file
	[data] = importfile(filename, 7, 8,['%*s%f%*[^\n]'],'\t');
	data   = cell2mat(data);
	if data(1)<0, error('"fs" must be positive'); end
	if data(2)<0, error('The variable "Duration" must be positive'); end

	fs   = data(1); % sampling frequency
	tmax = data(2); % duration of time series
	dt   = 1./fs; % time step, periodo
	N    = tmax/dt+1; % number of time step for a whole wind histories
	% build vector time and frequency
	t     = linspace(0,tmax,N);
	f0    = 1./tmax;
	if f0>fs,error('"fs" must be higher than 1/Duration');end
	f     = [f0:f0:fs/2];
	Nfreq = numel(f);
	df    = median(diff(f));

	%% WIND DATA
	[data] = importfile(filename, 12, 23,['%*s%f%*[^\n]'],'\t');
	data   = cell2mat(data);
	if any(data(1:3)  <  0), error('The standard deviations of wind velocity must be larger than 0'); end
	if any(data(4:6)  <= 0), error('The turbulence length scales must be positiv'); end
	if any(data(7:12) <  0), error('The decay coefficients must be positiv'); end
	
	sstdU = data(1);
	sstdV = data(2);
	sstdW = data(3);
	LLux  = data(4);
	LLvx  = data(5);
	LLwx  = data(6);
	Cuy   = data(7);
	Cuz   = data(8);
	Cvy   = data(9);
	Cvz   = data(10);
	Cwy   = data(11);
	Cwz   = data(12);



	[typeWind]  = importfile(filename, 29, 30,['%*s%s%*[^\n]'],'\t');
	type        = char(typeWind{1});
	windProfile = char(typeWind{2}); % power or log
	
	[data]    = importfile(filename, 32,36,['%*s%f%*[^\n]'],'\t');
	data      = cell2mat(data);
	Uref      = data(1);
	zr        = data(2);
	a         = data(3);
	u_star    = data(4);
	roughness = data(5);

	%% GRID GENERATION
	clear data
	[data] = importfile(filename, 40,45,['%*s%f%*[^\n]'],'\t');
	data   = cell2mat(data);
	if any(data(1:2)<0), error('Nyy and Nzz must be positiv'); end
	Nyy    = data(1);
	Nzz	   = data(2);
	Zmin   = data(3);
	Zmax   = data(4);
	Ymin   = data(5);
	Ymax   = data(6);

	if Zmin>Zmax,
		warning('Zmin > Zmax. Their values have been switched');
		dummy = Zmax;
		Zmax  = Zmin;
		Zmin  = dummy;
		clear dummy
	end

	if Ymin>Ymax,
		warning('Ymin > Ymax. Their values have been switched');
		dummy = Ymax;
		Ymax  = Ymin;
		Ymin  = dummy;
		clear dummy
	end

	% Check compatibility between grid and node number
	if and(Ymin == Ymax,Nyy>1),
		warning('Ymin = Ymax but Nyy > 1, Nyy is set to 1')
		Nyy = 1;
	end

	if and(Zmin == Zmax,Nzz>1),
		warning('Zmin = Zmax but Nzz > 1, Nzz is set to 1')
		Nzz = 1;
	end

% Create the grid
	clear nodes
	y       = linspace(Ymin,Ymax,Nyy);
	z       = linspace(Zmin,Zmax,Nzz);
	[Z,Y]   = meshgrid(z,y);
	nodes.Y = Y(:);
	nodes.Z = Z(:);



	% Wind profile
	if strcmpi(windProfile,'power'),
		U = Uref.*(nodes.Z./zr).^(a);
	elseif strcmpi(windProfile,'log'),
		k = 0.4; % Von karman constant
		U = u_star./k.*log(z./roughness);
	else
		error('wind profile selected is unknown\n')
	end
	nodes.U = U(:);
	% names affected at each nodes
	for ii = 1:Nyy*Nzz,
		nodes.name{ii} = strcat('N',num2str(ii));
	end

	stdU = zeros(Nzz,1);
	stdV = zeros(Nzz,1);
	stdW = zeros(Nzz,1);

	stdU(1) = sstdU;
	stdV(1) = sstdU;
	stdW(1) = sstdU;

	Lux = zeros(Nzz,1);
	Lvx = zeros(Nzz,1);
	Lwx = zeros(Nzz,1);

	Lux(1) = LLux;
	Lvx(1) = LLvx;
	Lwx(1) = LLwx;

	u_f    = zeros(Nzz,1);
	u_f(1) = u_star;
	
	v = 0.67+0.05*log(roughness);
	for i = 2:Nzz  
		u_f(i)  = (0.4*nodes.U(i))/(log(nodes.Z(i)/roughness)); 
		stdU(i) = (5.7*(u_f(i)^2))^(1/2); 
		stdV(i) = 0.75*stdU(i);
		stdW(i) = 0.5*stdU(i);
		Lux(i)  = 300*((nodes.Z(i)/200)^v);         
		Lvx(i)  = 0.3*Lux(i);
		Lwx(i)  = 0.2*Lux(i);
	end
	%% Input data for wind coherence
	dy    = zeros(Nyy*Nzz,Nyy*Nzz); % matrix distance along y
	dz    = zeros(Nyy*Nzz,Nyy*Nzz); % matrix distance along z
	meanU = zeros(Nyy*Nzz,Nyy*Nzz); % mean wind speed between two nodes
	for kk = 1:Nyy*Nzz,
		for mm = 1:Nyy*Nzz,
			dy(kk,mm)    = abs(nodes.Y(kk)-nodes.Y(mm));
			dz(kk,mm)    = abs(nodes.Z(kk)-nodes.Z(mm));
			meanU(kk,mm) =  0.5.*(nodes.U(kk)+nodes.U(mm));
		end
	end

	%% GENERATION OF WIND HISTORIES
	%preallocation
	clear ii jj A G dummy*
	dummySpeed = zeros(3*Nyy*Nzz,N);
	tic
	for jj = 1:Nfreq,
		%co-coherence
		DecY_u = dy.*Cuy*f(jj);
		DecY_v = dy.*Cvy*f(jj);
		DecY_w = dy.*Cwy*f(jj);
		
		DecZ_u = dz.*Cuz*f(jj);
		DecZ_v = dz.*Cvz*f(jj);
		DecZ_w = dz.*Cwz*f(jj);
		
		coh_u  = exp(-sqrt(DecY_u.^2+DecZ_u.^2)./meanU);
		coh_v  = exp(-sqrt(DecY_v.^2+DecZ_v.^2)./meanU);
		coh_w  = exp(-sqrt(DecY_w.^2+DecZ_w.^2)./meanU);
		
		% turbulence spectrum
		if strcmpi(type,'von karman'),
			Su = VonKarmanSpectrum(f(jj),nodes.U,stdU,Lux,'u');
			Sv = VonKarmanSpectrum(f(jj),nodes.U,stdV,Lvx,'v');
			Sw = VonKarmanSpectrum(f(jj),nodes.U,stdW,Lwx,'w');
		elseif strcmpi(type,'kaimal'),
			Su = KaimalSpectrum(f(jj),nodes.U,stdU,Lux,'u');
			Sv = KaimalSpectrum(f(jj),nodes.U,stdV,Lvx,'v');
			Sw = KaimalSpectrum(f(jj),nodes.U,stdW,Lwx,'w');
		else
			error(' spectrum type is unknown')
		end
		% spectral matrix with correlation between u and w
		Suu = sqrt(Su*Su').*coh_u;
		Svv = sqrt(Sv*Sv').*coh_v;
		Sww = sqrt(Sw*Sw').*coh_w;  
		S   = [Suu,zeros(size(Suu)),zeros(size(Suu));...
			zeros(size(Suu)),Svv,zeros(size(Suu));...
			zeros(size(Suu)),zeros(size(Suu)),Sww];
		phi = 2*pi.*rand(3*Nyy*Nzz,1);
		phi = repmat(phi,[1,N]);
		wt  = f(jj).*repmat(t,[3*Nyy*Nzz,1]);
		A   = cos(2.*pi.*wt+phi);
		G   = chol(S,'lower');% cholesky decomposition of Svv
		
		dummySpeed = dummySpeed + sqrt(2*df).*abs(G)*A;
	end
	toc
	% Output of the simulation
	u = dummySpeed(1:Nyy*Nzz,:);
	v = dummySpeed(Nyy*Nzz+1:2*Nyy*Nzz,:);
	w = dummySpeed(2*Nyy*Nzz+1:3*Nyy*Nzz,:);

	%% FUNCTIONS
    function [data] = importfile(filename, startRow, endRow,formatSpec,delimiter)
        % GOAL
        % extract data from txt files or excel files to store them into matrix
        
        %                               INPUT
        %  filename:
        %           type: string with extension.
        %           definition: name of the file whose information are extracted
        %  startRow:
        %           type: integer
        %           definition: first row of the file to be read
        %  endRow:
        %           type: integer
        %           definition: last row of the file to be read
        %  formatSpec:
        %           type: string
        %           definition: format of the data to be read.
        %           e.g. : formatSpec = ['%f%f%f%s%s%*[^\n]'];
        %  delimiter:
        %           type: string
        %           definition: symbol used to delimite the columns in the file
        %           e.g. : delimiter = '\t'; or delimiter = ',';
        
        %                               OUTPUT
        %  data:
        %           type: matrix [N1 x N2] where N1 and N2 are integers.
        %           N1 = endRow-startRow or total number of rows if enRow is not
        %           specified
        %           N2 = defined by formatSpec.
        %           definition: extracted data from the file, stored as a matrix
        % Copyright (C) Etienne Cheynet, 2015.
        % last modification: 27/01/2015 10:56
        % Initialize variables.
        if nargin<=2
            startRow = 1;
            endRow   = inf;
        end
        if endRow == [],
            endRow = inf;
        end
        %  open the text file.
        fileID = fopen(filename,'r');
        %Read columns of data according to format string.
        dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
        for block = 2:length(startRow)
            frewind(fileID);
            dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
            for col = 1:length(dataArray)
                dataArray{col} = [dataArray{col};dataArrayBlock{col}];
            end
        end
        % Close the text file.
        fclose(fileID);
        %% Create output variable
        dataArray = cellfun(@(x) num2cell(x), dataArray, 'UniformOutput', false);
        data = [dataArray{1:end}];
    end
    function [S] = VonKarmanSpectrum(f,V,stdV,L,component)
        % ---------------------------------------------
        % INPUT
        % f: float; frequency is [1 x 1]
        % V: float; Mean wind speed Normal to the deck is [1x1]
        % std_speed : float; std of speed is [1 x 1]
        % L =  float; turbulence length scales is [1x1]
        % stdV = ; float; std of wind velocity fluctuations [1x1]
        % component : string; is 'u' or 'w'
        % ---------------------------------------------
        % OUTPUT
        % Sv: float; [1x1] value of Spectrum for a given frequency
        % ---------------------------------------------
        % Von Karman coefficent
        coef = [-5/6, -11/6];
        % dimension of output
        S    = zeros(size(V)); % is [Nzz,1]
        %calculation of S/std^2
        n    = L.*V.^(-1).*f;
        if strcmpi(component,'u'),
            S =  V.^(-1).*4.*L.*stdV.^2.*(1+70.7.*n.^2).^(coef(1));
        elseif strcmpi(component,'v'),
            S =  V.^(-1).*4.*L.*stdV.^2.*(1+70.7.*4.*n.^2).^(coef(2)).*(1+188.4.*4*n.^2);
        elseif strcmpi(component,'w'),
            S =  V.^(-1).*4.*L.*stdV.^2.*(1+70.7.*4.*n.^2).^(coef(2)).*(1+188.4.*4*n.^2);
        else
            fprintf('error: component unknown \n\n')
            return
        end
    end
    function [S] = KaimalSpectrum(f,V,stdV,L,component)
        % ---------------------------------------------
        % INPUT
        % f: float; frequency is [1 x 1]
        % V: float; Mean wind speed Normal to the deck is [1x1]
        % std_speed : float; std of speed is [1 x 1]
        % Lturb =  float; turbulence length scales is [1x1]
        % stdV = ; float; std of wind velocity fluctuations [1x1]
        % component : string; is 'u' or 'w'
        % ---------------------------------------------
        % OUTPUT
        % Sv: float; [1x1] value of Spectrum for a given frequency
        % ---------------------------------------------
        % ---------------------------------------------
        %%
        n = f.*L./V;
        if strcmpi(component,'u'),
            S = stdV.^2./f.*(6.8.*n)./(1+1.5*6.8.*n).^(5/3);
        elseif strcmpi(component,'v'),
            S = stdV.^2./f.*(9.4.*n)./(1+1.5*9.4.*n).^(5/3);
        elseif strcmpi(component,'w'),
            S = stdV.^2./f.*(9.4.*n)./(1+1.5*9.4.*n).^(5/3);
        else
            fprintf(' spectrum type is unknown \n')
            return
        end
    end
end

