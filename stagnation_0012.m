clear all
NACA_0012;

%Create Angle of Attack
AoA = -16:.1:16;
AoA_rad = AoA.*(pi/180);
L = length(AoA_rad);

for i = 1:L
    [C_L_all_n(i,1),X_SEP_TOP(i,1)] = project3_code(AoA_rad(i));
end
X_SEP = X_SEP_TOP.*100;

L = length(X_SEP);
for i = 1:L
    if X_SEP(i,1)>20 && X_SEP(i+1,1)<20
       iSEP_bef_top = i;
       iSEP_aft_top = i+1;
       break;
    else
       continue;
    end
end

%Find theta at the two control points wrt TE bottom
aoa_sep_bef = AoA(iSEP_bef_top);
aoa_sep_aft = AoA(iSEP_aft_top);

%Linear Interpolate to Find Exact Angle of Bot Seprataion wrt TE
AOA_SEP_TOP = interp1([X_SEP(iSEP_bef_top,1);X_SEP(iSEP_aft_top,1)]...
            ,[aoa_sep_bef;aoa_sep_aft],20,'linear');

disp(AOA_SEP_TOP)
X_SEP(14) = [];
AoA(14) = [];

hold on
ylim([-10 100])
yline(20,'--r')
plot(AoA, X_SEP)
legend('Stall Angle Criteria','Separation Points','Location','Best')
title('NACA 0012 Separation Points vs Angle of Attack','Interpreter','latex')
xlabel('Angle of Attack ($\alpha$) [Degrees]','Interpreter','latex')
ylabel('$Sepration Point X Position [% of Chord]$','Interpreter','latex')