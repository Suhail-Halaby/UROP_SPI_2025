%% Suhail Halaby -- General Analysis Tool for RC Aircraft %% 
clc
clear 
close all
% Contents
    % 1 --> Vortex Panel Method Visualizer
    % 2 -->


%% Vortex Panel Method -- Analyzing 2D Airfoils

% Inputs 

    % Subdivisions (Panel Number)
    N = 200;
    % chord (m)
    c = 1;
    % freestream velocity (m/s) [X component, Y component]
    V_inf = 10;
    alpha = 20;
    % pressure conditions
    P_ambient = 101325;
    rho = 1.225;



    V_freestream = [V_inf*cosd(alpha),V_inf*sind(alpha)];



    % Airfoil Polar
    RawData = readmatrix("naca4415.dat");
    RawData(50,2) = 0;
    Xp =  c*RawData(2:end,1); 
    Yp =  c*RawData(2:end,2);
    clear RawData;

% Interpolating 

    % Upper Surface
    Yupr = Yp(Yp >= 0);
    Xupr = Xp(Yp >= 0);

    [Xupr,i_upr] = sort(Xupr);
    Xupr = [0;Xupr];
    Yupr = [0;Yupr(i_upr)];
    clear i_upr

    A = unique([Xupr,Yupr],"rows",'stable');
    clear Xupr Yupr
    Xupr = A(:,1);
    Yupr = A(:,2);
    clear A

    % Lower Surface
    Ylwr = Yp(Yp <= 0); 
    Xlwr = [Xp(Yp <= 0)];

    [Xlwr,i_lwr] = sort(Xlwr);
    Xlwr = [0;Xlwr];
    Ylwr = [0;Ylwr(i_lwr)];
    clear i_lwr

    B = unique([Xlwr,Ylwr],"rows",'stable');

    clear Xlwr Ylwr
    Xlwr = B(:,1);
    Ylwr = B(:,2);
    clear B




    % Cosine Spacing and Interpolation
    Xsp = c*(1 - cos(pi*(0:N)./N))/2;
    Yusp = interp1(Xupr,Yupr,Xsp,"spline");
    Ylsp = interp1(Xlwr,Ylwr,Xsp,"spline");

    % For Evaluating Nomal Vectors (type 2)
    Nfactor = 100;
    Nf = Nfactor*N;

    Xspc = c*(1 - cos(pi*(0:Nf)./Nf))/2;
    Yupc = interp1(Xupr,Yupr,Xspc,"spline");
    Ylwc = interp1(Xlwr,Ylwr,Xspc,"spline");
    
    clear Xlwr Ylwr Xupr Yupr Nf
    
    % Concatonnating Lists and Removing Redundant Points
    X = [Xsp';Xsp'];
    Y = [Yusp';Ylsp'];
    Xcurve = [Xspc';Xspc'];
    Ycurve = [Yupc';Ylwc'];
    

    % fix: need to loop from TE to TE, retain two panels at TE and LE, down
    % CCW - moves along top to bottom from the TE

    Xflip = flip(Xsp);
    Yflip = flip(Yusp);
    Xflipc = flip(Xspc);
    Yflipc = flip(Yupc);

    P = [Xflip(1:2:end)',Yflip(1:2:end)';Xsp(3:2:end)',Ylsp(3:2:end)'];
    C = [Xflip(2:2:end)',Yflip(2:2:end)';Xsp(2:2:end)',Ylsp(2:2:end)'];

    H = [Xflipc',Yflipc';Xspc',Ylwc'];


    
    Cpl1 = [Xflipc(2+Nfactor:2*Nfactor:end)',Yflipc(2+Nfactor:2*Nfactor:end)';Xspc(2+Nfactor:2*Nfactor:end)',Ylwc(2+Nfactor:2*Nfactor:end)'];
    Cmn1 = [Xflipc(Nfactor:2*Nfactor:end)',Yflipc(Nfactor:2*Nfactor:end)';Xspc(Nfactor:2*Nfactor:end)',Ylwc(Nfactor:2*Nfactor:end)'];

    % normal vectors
    dx = Cpl1(:,1)-Cmn1(:,1);
    normal = zeros(N,2);
    for k = 1:N
        if dx(k) > 0
           dy = Cpl1(k,2) - Cmn1(k,2); 
           normal(k,:) = [dy , -dx(k)]  / norm([dx(k),dy]) ;
        else
           dy = Cmn1(k,2) - Cpl1(k,2) ;
           normal(k,:) = [-dy , -dx(k)] / norm([dx(k),dy]);
        end
    end

    % there are now two more vorticies than there are collocation points
    % matrix A must be square, and the kutta condition may be implemented

    % Gamma --> N+2 vortex strengths --> {N+2} X {1} --> {N+1} X {1}
    % without LE vortex
    % A (inner) --> the impact of each of the N+2 vorticies on the N
    % collocation points --> {N} x {N+1}
    % the kutta condition can be appended onto this matrix 
    % and one can get rid of the additional LE votrex
    % A becomes --> {N+1} X {N+1}
    % RHS --> freestream and zero criculation stipulations --> {N+1} X {1}

    % Populate RHS

    RHS = zeros(N+1,1);
    for k = 1:N
        RHS(k) = dot(V_freestream,normal(k,:));
    end

    % Populate Inner A
    A = zeros(N+1);
    for i = 1:N % each row is the eqn corresponding to each collocation point
        for j = 1:N+1 % each col corresponds to each vortex's impact
            r = norm(C(i,:)-P(j,:));
            u = -( 1/(2*pi*r^2) )  *  (C(i,2)-P(j,2));
            v = ( 1/(2*pi*r^2) ) * (C(i,1)-P(j,1));
            A(i,j) = dot([u,v],normal(i,:));
        end
    end
    % Imposing Kutta Condition
    A(N+1,1) = 1;
    A(N+1,N+1) = 1;
    % Solve the System
    Gamma = pinv(A) * RHS;
    Cl = (2*sum(Gamma))/(c*norm(V_freestream))


    % Plotting the Velocity Feild
    reso = 0.01;
    SP = 5;
    M = 0.1*alpha+c;
    [Xview,Yview] = meshgrid(-c:reso:2*c,-M*c:reso:M*c);

    % Evaluating total velocity

    Velocity = zeros(size(Xview,1)*size(Xview,2),2);

    for i = 1:(size(Xview,1)*size(Xview,2))
        V_induced = zeros(1,2);
        for j = 1:N+1
            r = norm([Xview(i),Yview(i)]-P(j,:));
            V_induced(1) = V_induced(1) + ( Gamma(j)/(2*pi*r^2) ) * (Yview(i)-P(j,2));
            V_induced(2) = V_induced(2) - ( Gamma(j)/(2*pi*r^2) ) * (Xview(i)-P(j,1));
        end 
        Velocity(i,:) = V_freestream + V_induced;
    end

    Velocity(isnan(Velocity)) = 0;
    X_Vel = reshape(Velocity(:,1),size(Xview,1),size(Xview,2));
    Y_Vel = reshape(Velocity(:,2),size(Xview,1),size(Xview,2));







    
    subplot(3,1,1)
    hold on
    plot(P(:,1),P(:,2),'*r',MarkerEdgeColor='r',LineWidth=2,MarkerSize=6)
    plot(C(:,1),C(:,2),'*r',MarkerEdgeColor='b',LineWidth=2,MarkerSize=6)
    plot(H(:,1),H(:,2),'*r',MarkerEdgeColor='k',LineWidth=2,MarkerSize=2)
    plot(Cpl1(:,1),Cpl1(:,2),'*r',MarkerEdgeColor='g',LineWidth=2,MarkerSize=2)
    plot(Cmn1(:,1),Cmn1(:,2),'*r',MarkerEdgeColor='m',LineWidth=2,MarkerSize=2)
    streamline(Xview,Yview,X_Vel,Y_Vel,-c*ones(size(-M*c:reso*SP:M*c)),-M*c:reso*SP:M*c,Color='k')
    axis([-c 2*c -0.5*c 0.5*c])
    ax = gca;
    ax.DataAspectRatio = [1 1 1];

    subplot(3,1,2) 
    Levels = min(min((X_Vel.^2+Y_Vel.^2).^0.5)):0.5:max(max((X_Vel.^2+Y_Vel.^2).^0.5));
    contourf(Xview,Yview,(X_Vel.^2+Y_Vel.^2).^0.5,Levels)
    colorbarHandle = colorbar;
    ylabel(colorbarHandle, 'Velocity (m/s)');
    axis([-c 2*c -0.5*c 0.5*c])
    ax = gca;
    ax.DataAspectRatio = [1 1 1];


    subplot(3,1,3) 
    
    Total_Pressure = P_ambient + 0.5*rho*norm(V_freestream)^2;
    Static_Pressure = Total_Pressure - 0.5*rho*(X_Vel.^2+Y_Vel.^2);
    Levels = 0.99*Total_Pressure:20:max(max(Static_Pressure));



    contourf(Xview,Yview,Static_Pressure,Levels)
    colorbarHandle = colorbar;
    ylabel(colorbarHandle, 'Static Pressure (Pa)');
    axis([-c 2*c -0.5*c 0.5*c])
    ax = gca;
    ax.DataAspectRatio = [1 1 1];







