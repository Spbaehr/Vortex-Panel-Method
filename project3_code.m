function [C_L,X_SEP_TOP] = project3_code(alpha)

%Import geometry-------------------------------------------
%----------------------------------------------------------
Geometry = csvread("geometry.csv");

xbound = Geometry(:,1);
ybound = Geometry(:,2);
xcon = Geometry(:,3);
ycon = Geometry(:,4);
S = Geometry(:,5);

    %xbound = x coordinates of panel endpoints
    %ybound = y coordinates of panel endpoints
    %xcon = x coordinates of control points
    %ycon = y coordinates of control points
    %S = length of all panels
%End of Importing Geometry---------------------------------
%----------------------------------------------------------

%Numerical Pressure Coefficent Calculations----------------
%----------------------------------------------------------

N = length(Geometry) - 1;     %Get number of panels
THETA = zeros(N,1);           %Initialize theta values for panels
% alpha = 9*pi/180;

%Calculate Theta Values
for k = 1:N
    THETA(k) = atan2((ybound(k+1)-ybound(k)), (xbound(k+1)-xbound(k)));
end

%Calcul;ate Basic Trig Function Values
sine = sin(THETA);
cosine = cos(THETA);

%Caclulate RHS
RHS = sin(THETA - alpha);
RHS(N+1) = 0;

%Calulate Coefficient Matrices Cn1,Cn2,Ct1,Ct2-------------
%Initialize Coefficent Matrices
Cn1 = zeros(N,N);
Cn2 = zeros(N,N);
Ct1 = zeros(N,N);
Ct2 = zeros(N,N);

for i = 1:N
    for j = 1:N
        if i == j
           Cn1(i,j) = -1;
           Cn2(i,j) = 1;
           Ct1(i,j) = 0.5*pi;
           Ct2(i,j) = 0.5*pi; 
        else
           A = -((xcon(i)-xbound(j)).*cosine(j))...
               -((ycon(i)-ybound(j)).*sine(j));
           B = (xcon(i)-xbound(j)).^2 + (ycon(i)-ybound(j)).^2;
           C = sin(THETA(i)-THETA(j));
           D = cos(THETA(i)-THETA(j));
           E = ((xcon(i)-xbound(j)).*sine(j))...
               -((ycon(i)-ybound(j)).*cosine(j));
           F = log(1 + (S(j).*(S(j)+2.*A))./B);
           G = atan2(E.*S(j), B+A.*S(j));
           P = (xcon(i)-xbound(j)).*sin(THETA(i)-2.*THETA(j))...
                +(ycon(i)-ybound(j)).*cos(THETA(i)-2.*THETA(j));
           Q = (xcon(i)-xbound(j)).*cos(THETA(i)-2.*THETA(j))...
               -(ycon(i)-ybound(j)).*sin(THETA(i)-2.*THETA(j));
           Cn2(i,j) = D + (0.5.*Q.*F)./S(j) - ((A.*C+D.*E).*G)./S(j); 
           Cn1(i,j) = 0.5.*D.*F + C.*G - Cn2(i,j);
           Ct2(i,j) = C + (0.5.*P.*F)./S(j) + ((A.*D-C.*E).*G)./S(j);
           Ct1(i,j) = 0.5.*C.*F - D.*G - Ct2(i,j);
        end
    end
end
%----------------------------------------------------------

%Calulate Influence Coefficients Matrices An,At------------
An = zeros(N+1,N+1);    %Initialize An matrix
At = zeros(N,N+1);      %Initialize At matrix

for a = 1:N
    An(a,1) = Cn1(a,1);
    An(a,N+1) = Cn2(a,N);
    At(a,1) = Ct1(a,1);
    At(a,N+1) = Ct2(a,N); 
    for b = 2:N
        An(a,b) = Cn1(a,b) + Cn2(a,b-1);
        At(a,b) = Ct1(a,b) + Ct2(a,b-1);
    end
end

An(N+1,1) = 1; 
An(N+1,N+1) = 1;
for c = 2:N
   An(N+1,c) = 0; 
end
%----------------------------------------------------------

%Calculate Gamma (Vortex Panel Strength)
gamma = An\RHS;

%Calculate Velocity/Pressure Coefficients at Control Points
Ut = zeros(N,1);      %Initialize velocity coef. matrix
C_p = zeros(N,1);     %Initialize Pressure coef. matrix

for d = 1:N
    Ut(d) = cos(THETA(d) - alpha);
    for e = 1:N+1
       Ut(d) = Ut(d) + At(d,e).*gamma(e);
       C_p(d) = 1 - Ut(d).^2;
    end
end

%End of Numerical Pressure Coefficient Calculations--------
%----------------------------------------------------------

%Plot Pressure Coefficients--------------------------------
%----------------------------------------------------------
% fig1 = figure;
% c = 1;
% x_c_L = xcon(1:end/2)/c;
% x_c_U = xcon(end/2:end-1)/c;
% %Align Cp values with correct theta values
% Cp_L = C_p(1:end/2);
% Cp_U = C_p(end/2+1:end);
% 
% hold on
% xlim([0 1.1])
% ylim([-17 1.5])
% plot(x_c_L,Cp_L,'-o')
% plot(x_c_U,Cp_U,'-o')
% title('NACA 0012 $C_{p}$ for N=64 Panels and AoA=$-16^{\circ}$','Interpreter','latex')
% legend('$C_{p}$ Lower Surface',' $C_{p}$ Upper Surface', 'Location','best','Interpreter','latex')
% xlabel('X/C')
% ylabel('$C_{p}$','Interpreter','latex')
% hold off

%End of Plot Pressure Coefficients-------------------------
%----------------------------------------------------------

%Set Constant Parameters
P_inf = 101300;      %Pressure [Pa]
rho = 1.225;         %Density [kg/m^3]
U_inf = 25;         %Velcoity [m/s]

%Convert Cp values to actual pressures P
for j = 1:length(C_p)
    P(j) = (0.5.*C_p(j).*rho*U_inf.^2) + P_inf;
    dP(j) = (0.5.*C_p(j).*rho*U_inf.^2);
end

%Calculate theta wrt positive x and second point
Theta = atan2((ybound(2:end)-ybound(1:end-1)), ...
              (xbound(2:end)-xbound(1:end-1)));
Theta_d = Theta.*(180/pi);

%Calculate normal vectors
nX = -sin(Theta);
nY = cos(Theta);

% plot(xcon,ycon,'ko')
% quiver(xcon(1:end-1),ycon(1:end-1),nX,nY,'color','bl')

%Calculate forces per unit span
for j = 1:length(P)
    fX(j) = -P(j).*nX(j).*S(j);
    fY(j) = -P(j).*nY(j).*S(j);
end

%Calculate net forces
FX = sum(fX);
FY = sum(fY);

%Caclulate forces in lift and drag directions
T = [cos(alpha),sin(alpha);-sin(alpha),cos(alpha)];
F = [FX;FY];
F_DL = T*F;

C_L = F_DL(2)/(0.5*rho*(U_inf^2));

%Separation Point Calculations-----------------------------
%----------------------------------------------------------
%Find Leading Edge
indexLE = find(Ut(1:end-1)<0 & Ut(2:end)>0);

%Flip Ut so in order from LE and Split into two sections
Ut_new_bot = flip(Ut(1:indexLE,1));    
Ut_new_top = Ut(indexLE+1:end,1);

%Create distance traveled vector----------------------------
L_bot = length(S(1:indexLE,1));
L_top = length(S(indexLE+1:end-1,1));
x_s_bot = zeros(L_bot,1);
x_s_top = zeros(L_top,1);
S(1:indexLE,1) = flip(S(1:indexLE,1));
s_dis_bot = S(1,1)/2;           %Initial distance starting at first control point
s_dis_top = S(indexLE+1,1)/2;   %Initial distance starting at first control point

%For Bottom
for p = 1:L_bot
    x_s_bot(p,1) = s_dis_bot;
    s_dis_bot = s_dis_bot + S(p+1,1);
end

%For Top
for p = indexLE+1:L_top+indexLE
    x_s_top(p-L_bot,1) = s_dis_top;
    s_dis_top = s_dis_top + S(p+1,1);
end
%----------------------------------------------------------

%Integration of Ut
int_Ut_bot = cumtrapz(x_s_bot,Ut_new_bot.^5);    %Integrate along x_s
int_Ut_top = cumtrapz(x_s_top,Ut_new_top.^5);    %Integrate along x_s

%Differentiation of Ut
dy_dx_bot = zeros(L_bot,1);
dy_dx_top = zeros(L_top,1);
for i = 2:L_bot-1
    dy_dx_bot(i-1,1) = (Ut_new_bot(i+1,1)...
                       -Ut_new_bot(i-1,1))/(2*(x_s_bot(i,1)-x_s_bot(i-1,1))); %CD formula 
end

for i = 2:L_top-1
    dy_dx_top(i,1) = (Ut_new_top(i+1,1)...
                     -Ut_new_top(i-1,1))/(2*(x_s_top(i,1)-x_s_top(i-1,1))); %CD formula
end

%Evaluating K
K_bot = (0.45./(Ut_new_bot(2:end-1,1).^6)).*dy_dx_bot(2:end-1,1)...
        .*int_Ut_bot(2:end-1,1);
K_top = (0.45./(Ut_new_top(2:end-1,1).^6)).*dy_dx_top(2:end-1,1)...
        .*int_Ut_top(2:end-1,1);
    
%Find Theta of Separation Points Bottom-----------------------
L = length(K_top);
for i = 1:L
    if K_top(i,1)>-0.09 && K_top(i+1,1)<-0.09
       iSEP_bef_top = i;
       iSEP_aft_top = i+1;
       break;
    else
       continue;
    end
end

%Find theta at the two control points wrt TE bottom
x_sep_top_bef = xbound(iSEP_bef_top+indexLE+1,1);
x_sep_top_aft = xbound(iSEP_aft_top+indexLE+1,1);

%Linear Interpolate to Find Exact Angle of Bot Seprataion wrt TE
X_SEP_TOP = interp1([K_top(iSEP_bef_top,1);K_top(iSEP_aft_top,1)]...
            ,[x_sep_top_bef;x_sep_top_aft],-0.09,'linear');
%----------------------------------------------------------

%End of Separation Point Calculations----------------------
%----------------------------------------------------------
end