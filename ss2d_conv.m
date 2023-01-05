function [xc,yc,xf,yf,A,b,U] = ss2d_conv(inputFile)

%% Import problem file
[vars, vals] = importInputs();

%% Extract input values

% Geometry
t = str2double(vals{contains(vars,'thickness')});       % Get domian thickness
Lx = str2double(vals{contains(vars,'Lx')});     % Get domian length in X
Ly = str2double(vals{contains(vars,'Ly')});     % Get domian length in Y
dx = vals(contains(vars,'dx'),:);                 % Resolution in X
dy = vals(contains(vars,'dy'),:);                 % Resolution in Y

% Materials
materials = vals(contains(vars,'material')',:);

% Boundary Conditions
bcx_min = vals(contains(vars,'xmin')',:);     % BC in min X
bcx_max = vals(contains(vars,'xmax')',:);     % BC in max X
bcy_min = vals(contains(vars,'ymin')',:);     % BC in max Y
bcy_max = vals(contains(vars,'ymax')',:);     % BC in max Y

% Generation terms
sources = vals(contains(vars,'source')',:);

% Convection control
convection = vals(contains(vars,'convection')',:);     % convection

%% Generate mesh

[xc, xf] = gen1DmeshVals(Lx,dx,'x');    % Genrate x cell and face coords
[yc, yf] = gen1DmeshVals(Ly,dy,'y');    % Genrate x cell and face coords
coords_cells = genCoords2dArray(xc,yc);  % Array of cell coordinate pairs
coords_faces = genCoords2dArray(xf,yf);  % Array of face coordinate pairs
N = size(coords_cells,1);               % Amount of cells

%% Formulate problem: Generate A Matrix and B Vector

[A,b] = genModel();

%% Solve model

U = A\b;

%% Visualizations

figure('units','normalized','outerposition',[0 0 1 1]);
[XC,YC,XF,YF] = plotMesh();
figure('units','normalized','outerposition',[0 0 1 1]);
plotResults();
figure('units','normalized','outerposition',[0 0 1 1]);
plotResultsAndMesh();

%% Functions

    function [vars, vals] = importInputs()
        fid = fopen(inputFile); i = 1;
        vals = cell(1,20);
        while ~feof(fid)

            % Read string and split by commas
            strelems = strsplit(strip(fgets(fid)),',');
            if ~isempty(strelems{1})
                if length(strelems) > 1
                    vars(i,1) = strelems(1);
                    vals(i,1:length(strelems(2:end))) = strelems(2:end);
                    i = i + 1;
                end
            end

        end
        fclose(fid);
    end

    function [dc,df] = gen1DmeshVals(L,def,flag)
        res = def{1};
        if strcmp(res,'linear')

            % Check if generating x, y or z coordinates
            if ~isempty(flag)
                if strcmp(flag,'x')
                    seed = str2double(def{2});
                    seedf = str2double(def{3});
                elseif strcmp(flag,'y')
                    seed = str2double(def{2});
                    seedf = str2double(def{3});
                else
                    seed = str2double(def{2});
                    seedf = str2double(def{3});
                end
            end

            % Initialize dc and df
            df = 0;

            % Generate face coordinates
            i = 2;
            while df(end) < L
                if (i-1)*seed >= seedf
                    df(i) = df(i-1)+seedf; % Face coordinates
                    dc(i-1) = 0.5*(df(i-1)+df(i)); % Cell coordinates
                    i = i + 1;
                else
                    df(i) = df(i-1)+(i-1)*seed; % Face coordinates
                    dc(i-1) = 0.5*(df(i-1)+df(i)); % Cell coordinates
                    i = i + 1;
                end
            end

            % Correct the last coordinate
            df(end) = L;
            dc(end) = 0.5*(df(end) + df(end-1));

        elseif strcmp(res,'quadratic')

        elseif strcmp(res,'expoential')

        elseif strcmp(res,'usr_def')

        else
            delta = str2double(res);
            d = 0:delta/2:L;
            dc = d(2:2:end-1);
            df = d(1:2:end);
        end

    end

    function [XC,YC,XF,YF] = plotMesh()
        [XF,YF] = meshgrid(xf,yf);
        [XC,YC] = meshgrid(xc,yc);
        mesh(XF,YF,zeros(size(XF))); view ([0 0 90]); hold on;
        plot3(XC,YC,zeros(size(XC)),'.b');
        xlabel('X'); ylabel('Y'); zlabel('Z'); axis tight;
        title('Mesh'); legend('Cell Faces','Cell Centers');
    end

    function [] = plotResults()
        contourf(XC', YC', reshape(U,length(xc),length(yc)));
        colorbar; title('Results'); xlabel('X'); ylabel('Y');
        colormap(jet);
    end

    function [XC,YC,XF,YF] = plotResultsAndMesh()
        [XF,YF] = meshgrid(xf,yf);
        [XC,YC] = meshgrid(xc,yc);
        mesh(XF,YF,zeros(size(XF))); view ([0 0 90]); hold on;
        plot3(XC,YC,zeros(size(XC)),'.b','MarkerSize',20);
        xlabel('X'); ylabel('Y'); zlabel('Z'); axis tight;
        legend('Cell Faces','Cell Centers');
        hold on; surf(XC', YC', reshape(U,length(xc),length(yc)),...
            'FaceAlpha',0.5,'EdgeColor','None');
        colorbar; title('Results'); xlabel('X'); xlabel('Y');
        title('Results over the Mesh');
        colormap(jet);
    end

    function coords = genCoords3dArray(x,y,z)
        count = 1;
        for i = 1:length(x)
            for j = 1:length(y)
                for k = 1:length(z)
                    coords(count,:) = [x(i) y(j) z(k)];
                    count = count + 1;
                end
            end
        end
    end

    function coords = genCoords2dArray(x,y)
        count = 1;
        for i = 1:length(x)
            for j = 1:length(y)
                coords(count,:) = [x(i) y(j)];
                count = count + 1;
            end
        end
    end

    function [A,b] = genModel()

        counter = 1; b = zeros(length(xc)*length(yc),1);
        convScheme = convection{1};
        vx = str2double(convection{2});
        vy = str2double(convection{3});
        for j = 1:length(yc)
            for i = 1:length(xc)

                % Compute face deltas
                deltaXf = xf(i+1) - xf(i);
                deltaYf = yf(j+1) - yf(j);

                % Compute heat source
                cellArea = deltaXf*deltaYf;
                cellVol = cellArea*t;

                % Get sources
                for s = 1:size(sources,1)

                    % Get source data values
                    source = sources(s,:);

                    % Get the domain coordinates of the source
                    xs = str2double(source{2});
                    xe = str2double(source{3});
                    ys = str2double(source{4});
                    ye = str2double(source{5});

                    % Solve the source if the current cell is withing the source
                    % domain.
                    if xf(i) >= xs && xf(i+1) <= xe && yf(j) >= ys && yf(j+1) <= ye
                        func = eval(['@(x,y) ' source{1}]);
                        b(counter,1) = b(counter,1) +...
                            cellVol*func(xc(i),yc(j));
                    end
                end

                % Get material properties
                for m = 1:size(materials,1)
                    material = materials(m,:);

                    % Get the domain coordinates of the material
                    xs = str2double(material{7});
                    xe = str2double(material{8});
                    ys = str2double(material{9});
                    ye = str2double(material{10});

                    % Get the density, spec heat and viscosity
                    rho = str2double(material{1});
                    cp = str2double(material{2});
                    visc = str2double(material{3});

                    if xf(i) >= xs && xf(i+1) <= xe && yf(j) >= ys && yf(j+1) <= ye

                        % Compute diffusion coefficient in X
                        func = eval(['@(x,y) ' material{4}]);
                        kx = func(xc(i),yc(j));

                        % Compute diffusion coefficient in Y
                        func = eval(['@(x,y) ' material{5}]);
                        ky = func(xc(i),yc(j));

                        % Compute diffusion coefficient in Z
                        %                         func = eval(['@(x,y) ' material{6}]);
                        %                         kz = func(xc(i),yc(j));
                    end

                end

                % Compute flows
                Qx =  vx*rho*deltaYf*t;
                Qy =  vy*rho*deltaXf*t;
                Qxp = Qx;
                Qxn = Qx;
                Qyp = Qy;
                Qyn = Qy;

                % If statements to check if at domain boundary in X
                if i == 1

                    % Interpolation Terms
                    sxp = (xf(i+1) - xc(i))/(xc(i+1) - xc(i));
                    sxm = (xf(i) - xc(i))/(xf(i) - xc(i));

                    % Compute coefficients
                    deltaXp = xc(i+1) - xc(i);
                    deltaXn = xc(i) - xf(i);
                    alphap = t*kx*deltaYf/deltaXp;

                    if strcmp(bcx_min{2},'Drichlet')
                        alphan = 2*t*kx*deltaYf/deltaXn;
                    elseif strcmp(bcx_min{2},'Neumann')
                        alphan = 2*t*kx*deltaYf;
                    else
                        fprintf('Incorrect BC definition in "xmin".\n');
                        fprintf('BC options are: "Drichlet" or "Neumann"\n');
                        fprintf('Exiting.\n'); return;
                    end

                    % Generate BC function
                    func = eval(['@(x,y) ' bcx_min{1}]);

                    % Compute Peclet NUmbers
                    Pexp = Qxp/alphap;
                    Pexn = Qxn/alphan;

                    % Populate A with coefficient
                    if strcmp(convScheme,'cd')
                        A(counter,counter+1) = Qxp*sxp - alphap;
                        b(counter,1) = b(counter,1) + func(xc(i),yc(j))*(alphan+Qxn*sxm);
                    elseif strcmp(convScheme,'upwind')
                        A(counter,counter+1) = max(0,-Qxp) - alphap;
                        b(counter,1) = b(counter,1) + func(xc(i),yc(j))*(alphan+max(0,Qxn));
                    elseif strcmp(convScheme,'hybrid')
                    end

                elseif i == length(xc)

                    % Interpolation Terms
                    sxp = (xf(i+1) - xc(i))/(xf(i+1) - xc(i));
                    sxm = (xf(i) - xc(i))/(xc(i-1) - xc(i));

                    % Compute coefficients
                    deltaXp = xf(end) - xc(i);
                    deltaXn = xc(i) - xc(i-1);
                    alphan = t*kx*deltaYf/deltaXn;

                    if strcmp(bcx_max{2},'Drichlet')
                        alphap = 2*t*kx*deltaYf/deltaXp;
                    elseif strcmp(bcx_max{2},'Neumann')
                        alphap = 2*t*kx*deltaYf;
                    else
                        fprintf('Incorrect BC definition in "xmax".\n');
                        fprintf('Exiting.\n'); return;
                    end

                    % Generate BC function
                    func = eval(['@(x,y) ' bcx_max{1}]);

                    % Compute Peclet NUmbers
                    Pexp = Qxp/alphap;
                    Pexn = Qxn/alphan;

                    % Populate A with coefficient
                    if strcmp(convScheme,'cd')
                        A(counter,counter-1) = -(alphan + Qxn*sxm);
                        b(counter,1) = b(counter,1) + func(xc(i),yc(j))*(alphap-Qxp*sxp);
                    elseif strcmp(convScheme,'upwind')
                        A(counter,counter-1) = -(alphan + max(0,Qxn));
                        b(counter,1) = b(counter,1) + func(xc(i),yc(j))*(alphap-max(0,-Qxp));
                    elseif strcmp(convScheme,'hybrid')
                    end

                else

                    % Interpolation Terms
                    sxp = (xf(i+1) - xc(i))/(xc(i+1) - xc(i));
                    sxm = (xf(i) - xc(i))/(xc(i-1) - xc(i));

                    deltaXp = xc(i+1) - xc(i);
                    deltaXn = xc(i) - xc(i-1);
                    alphap = t*kx*deltaYf/deltaXp;
                    alphan = t*kx*deltaYf/deltaXn;

                    % Compute Peclet NUmbers
                    Pexp = Qxp/alphap;
                    Pexn = Qxn/alphan;

                    % Compute coefficients
                    if strcmp(convScheme,'cd')
                        A(counter,counter+1) = Qxp*sxp - alphap;
                        A(counter,counter-1) = -(alphan + Qxn*sxm);
                    elseif strcmp(convScheme,'upwind')
                        A(counter,counter+1) = max(0,-Qxp) - alphap;
                        A(counter,counter-1) = -(alphan + max(0,Qxn));
                    elseif strcmp(convScheme,'hybrid')
                    end

                end

                % If statements to check if at domain boundary in Y
                if j == 1

                    % Interpolation Terms
                    syp = (yf(j+1) - yc(j))/(yc(j+1) - yc(j));
                    sym = (yf(j) - yc(j))/(yf(j) - yc(j));
                    %sym = 0.5;

                    % Compute coefficients
                    deltaYp = yc(j+1) - yc(j);
                    deltaYn = yc(j) - yf(j);
                    betap = t*ky*deltaXf/deltaYp;

                    if strcmp(bcy_min{2},'Drichlet')
                        betan = 2*t*ky*deltaXf/deltaYn;
                    elseif strcmp(bcy_min{2},'Neumann')
                        betan = 2*t*ky*deltaXf;
                    else
                        fprintf('Incorrect BC definition in "ymin".\n');
                        fprintf('Exiting.\n'); return;
                    end

                    % Generate BC function
                    func = eval(['@(x,y) ' bcy_min{1}]);

                    % Compute Peclet NUmbers
                    Peyp = Qyp/betap;
                    Peyn = Qyn/betan;

                    % Populate A with coefficient
                    if strcmp(convScheme,'cd')
                        A(counter,counter+length(xc)) = Qyp*syp - betap;
                        b(counter,1) = b(counter,1) + str2double(bcy_min{1})*(betan+Qyn*sym);
                    elseif strcmp(convScheme,'upwind')
                        A(counter,counter+length(xc)) = max(0,-Qyp) - betap;
                        b(counter,1) = b(counter,1) + func(xc(i),yc(j))*(betan+max(0,Qyn));
                    elseif strcmp(convScheme,'hybrid')
                    end

                elseif j == length(yc)

                    % Interpolation Terms
                    syp = (yf(j+1) - yc(j))/(yf(j+1) - yc(j));
                    %syp = 0.5;
                    sym = (yf(j) - yc(j))/(yc(j-1) - yc(j));

                    % Compute coefficients
                    deltaYp = yf(end) - yc(j);
                    deltaYn = yc(j) - yc(j-1);
                    betan = t*ky*deltaXf/deltaYn;

                    if strcmp(bcy_max{2},'Drichlet')
                        betap = 2*t*ky*deltaXf/deltaYp;
                    elseif strcmp(bcy_max{2},'Neumann')
                        betap = 2*t*ky*deltaXf;
                    else
                        fprintf('Incorrect BC definition in "ymax".\n');
                        fprintf('Exiting.\n'); return;
                    end

                    % Generate BC function
                    func = eval(['@(x,y) ' bcy_max{1}]);

                    % Compute Peclet NUmbers
                    Peyp = Qyp/betap;
                    Peyn = Qyn/betan;

                    % Coefficients
                    if strcmp(convScheme,'cd')
                        A(counter,counter-length(xc)) = -(betan + Qyn*sym);
                        b(counter,1) = b(counter,1) + str2double(bcy_max{1})*(betap-Qyp*syp);
                    elseif strcmp(convScheme,'upwind')
                        A(counter,counter-1) = -(betan + max(0,Qyn));
                        b(counter,1) = b(counter,1) + func(xc(i),yc(j))*(betap-max(0,-Qyp));
                    elseif strcmp(convScheme,'hybrid')
                    end

                else

                    % Interpolation Terms
                    syp = (yf(j+1) - yc(j))/(yc(j+1) - yc(j));
                    sym = (yf(j) - yc(j))/(yc(j-1) - yc(j));

                    % Compute coefficients
                    deltaYp = yc(j+1) - yc(j);
                    deltaYn = yc(j) - yc(j-1);
                    betap = t*ky*deltaXf/deltaYp;
                    betan = t*ky*deltaXf/deltaYn;

                    % Compute Peclet NUmbers
                    Peyp = Qyp/betap;
                    Peyn = Qyn/betan;

                    % Coeficients
                    if strcmp(convScheme,'cd')
                        A(counter,counter+length(xc)) = Qyp*syp - betap;
                        A(counter,counter-length(xc)) = -(betan + Qyn*sym);
                    elseif strcmp(convScheme,'upwind')
                        A(counter,counter+length(xc)) = max(0,-Qyp) - betap;
                        A(counter,counter-length(xc)) = -(betan + max(0,Qyn));
                    elseif strcmp(convScheme,'hybrid')
                    end
                end

                % Populate cell (i,j) coefficient
                if strcmp(convScheme,'cd')
                    A(counter,counter) = (alphap + alphan + betap + betan +...
                        Qxp*(1 - sxp) - Qxn*(1 - sxm) + Qyp*(1 - syp) - Qyn*(1 - sym));
                elseif strcmp(convScheme,'upwind')
                    A(counter,counter) = (alphap + alphan + betap + betan +...
                        max(0,Qxp) - max(0,-Qxn) + max(0,Qyp) - max(0,-Qyn));
                elseif strcmp(convScheme,'hybrid')
                end

                % Increase counter
                counter = counter + 1;
            end
        end
    end

end