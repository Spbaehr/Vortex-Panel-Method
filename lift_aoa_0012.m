clear all
NACA_0012;
%Create Angle of Attack
AoA = -16:.1:16;
AoA_rad = AoA.*(pi/180);

L = length(AoA_rad);
C_L_all = zeros(L,1);
for i = 1:L
    [C_L_all_n(i,1),X_SEP_TOP(i,1)] = project3_code(AoA_rad(i));
end

for i = 1:L
    C_L_all_a(i,1) = 2*pi*AoA_rad(i);
end

hold on
ylim([-2 2.5])
plot(AoA, C_L_all_n)
plot(AoA, C_L_all_a)
legend('Numerical Model','Analytical Model','Location','Best')
title('NACA 0012 $C_{L}$ vs Anghle of Attack','Interpreter','latex')
xlabel('Angle of Attack ($\alpha$) [Degrees]','Interpreter','latex')
ylabel('$C_{L}$','Interpreter','latex')

