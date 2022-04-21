clear all

n = 64;    %Number of Panels

%Create AIRFOIL--------------
theta(:,1) = linspace(-30,330,n+1);
bound =  [1.000000 -0.000000
          0.997638  0.000667
          0.990573  0.002647
          0.978864  0.005880
          0.962609  0.010271
          0.941947  0.015691
          0.917057  0.021985
          0.888157  0.028979
          0.855503  0.036484
          0.819389  0.044304
          0.780144  0.052241
          0.738132  0.060097
          0.693744  0.067680
          0.647402  0.074804
          0.599548  0.081294
          0.550646  0.086988
          0.501174  0.091737
          0.451623  0.095412
          0.402486  0.097904
          0.353518  0.098838
          0.305921  0.097840
          0.260260  0.094970
          0.217037  0.090347
          0.176728  0.084145
          0.139770  0.076589
          0.106553  0.067939
          0.077413  0.058482
          0.052633  0.048514
          0.032437  0.038325
          0.016990  0.028179
          0.006405  0.018304
          0.000739  0.008874
          0.000000  0.000000
          0.004076 -0.007913
          0.012810 -0.014507
          0.026069 -0.019799
          0.043684 -0.023825
          0.065446 -0.026642
          0.091117 -0.028326
          0.120437 -0.028982
          0.153123 -0.028733
          0.188878 -0.027732
          0.227393 -0.026150
          0.268344 -0.024177
          0.311395 -0.022012
          0.356197 -0.019857
          0.402423 -0.017905
          0.450360 -0.015990
          0.498826 -0.013960
          0.547371 -0.011922
          0.595542 -0.009966
          0.642883 -0.008158
          0.688939 -0.006542
          0.733265 -0.005140
          0.775426 -0.003957
          0.815005 -0.002983
          0.851604 -0.002197
          0.884854 -0.001576
          0.914412 -0.001092
          0.939974 -0.000721
          0.961271 -0.000443
          0.978077 -0.000242
          0.990212 -0.000105
          0.997546 -0.000026
          1.000000  0.000000];
xbound = bound(:,1);
ybound = flip(bound(:,2));

%Calculate the panel lengths--------------------------
S = zeros(length(xbound),1);

for i = 1:length(xbound)-1
    S(i,1) = sqrt((xbound(i+1,1)-xbound(i,1))^2 +...
                  (ybound(i+1,1)-ybound(i,1))^2);
end

%Calculate Control Points----------------------------
xcon = zeros(length(xbound),1);
ycon = zeros(length(xbound),1);

for i = 1:length(xbound)-1
    xcon(i,1) = (xbound(i+1,1) + xbound(i,1))/2;
    ycon(i,1) = (ybound(i+1,1) + ybound(i,1))/2;
end

%Congregate all values for export
Geometry = [xbound,ybound,xcon,ycon,S];

%Export Data as CSV File
writematrix(Geometry,'geometry.csv')