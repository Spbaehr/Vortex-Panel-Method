clear all
NACA_4412;
%Create Angle of Attack
AoA = -16:.1:16;
AoA_rad = AoA.*(pi/180);

L = length(AoA_rad);
C_L_all = zeros(L,1);
for i = 1:L
    C_L_all_n(i,1) = project3_code(AoA_rad(i));
end

for i = 1:L
    C_L_all_a(i,1) = 2*pi*(AoA_rad(i)- (-4.5*pi/180));
end

hold on
ylim([-2 2.5])
plot(AoA, C_L_all_n)
plot(AoA, C_L_all_a)
legend('Numerical Model','Analytical Model','Location','Best')
title('NACA 4412 $C_{L}$ vs Anghle of Attack','Interpreter','latex')
xlabel('Angle of Attack ($\alpha$) [Degrees]','Interpreter','latex')
ylabel('$C_{L}$','Interpreter','latex')