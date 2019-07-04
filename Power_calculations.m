%Akhil Kandhari, axk751@case.edu
%Review article math for design condition based on geometrical properties and wave properties

%Assuming all segments have same mass
%Assuming all segments can move at constant maximum linear speed

%Variables used
n = 6; %Total number of segments
eps_dot = 1; %Rate of change of length
g = 981; %Accelaration due to gravity in cm/s/s
Ro = 23; %Outer radius of segment in cm
Ri = 12; %Inner radius of segment in cm
L = 20; %Segment length in cm
nu = 1.3; %Poisson's ratio of structure

%Geometrical properties of structure:
I = (pi/4)*(Ro^4-Ri^4); %Second area moment of segment
V = pi*(Ro^2-Ri^2)*L; %Volume of a segment

geometry_factor = (g*V)/(2*nu*L^2);

%%
w = 1; %Number of waves occuring at one point
m = 1; %Pairs of maoving segments throughout the body
b = 0; %Number of bridged segments in between moving pairs

%Speed of a peristaltic robot:
max_theoretical_velocity = w*m*(m+b)*eps_dot/n;

%Calculation of squishing factor and squishing geometry
Ws = (5*n^2)/(4*(n-w*(2*m+b))*(w*m*(m+b)));
Gs = (Ro^2)/(Ro^3-Ri^3);

%Calculation of bending factor and bending geometry
Wb = ((2*m+b)*n)/(2*pi*w*m*(m+b));
Gb = (L^4)/(Ro*(Ro^4-Ri^4));

desired_velocity_2x1 = 0:.001:max_theoretical_velocity; %In cm/s
%Required power to weight ratio for given waves and geometry of structure and desired velocity
for i = 1:length(desired_velocity_2x1)
    power_to_weight_2x1(i) = (Ws*Gs+Wb*Gb)*geometry_factor*desired_velocity_2x1(i);
end

%%
w = 2; %Number of waves occuring at one point
m = 1; %Pairs of maoving segments throughout the body
b = 0; %Number of bridged segments in between moving pairs

%Speed of a peristaltic robot:
max_theoretical_velocity = w*m*(m+b)*eps_dot/n;

%Calculation of squishing factor and squishing geometry
Ws = (5*n^2)/(4*(n-w*(2*m+b))*(w*m*(m+b)));
Gs = (Ro^2)/(Ro^3-Ri^3);

%Calculation of bending factor and bending geometry
Wb = ((2*m+b)*n)/(2*pi*w*m*(m+b));
Gb = (L^4)/(Ro*(Ro^4-Ri^4));

desired_velocity_2x2 = 0:.001:max_theoretical_velocity; %In cm/s
%Required power to weight ratio for given waves and geometry of structure and desired velocity
for i = 1:length(desired_velocity_2x2)
    power_to_weight_2x2(i) = (Ws*Gs+Wb*Gb)*geometry_factor*desired_velocity_2x2(i);
end

%%
w = 1; %Number of waves occuring at one point
m = 1; %Pairs of maoving segments throughout the body
b = 1; %Number of bridged segments in between moving pairs
%Speed of a peristaltic robot:
max_theoretical_velocity = w*m*(m+b)*eps_dot/n;

%Calculation of squishing factor and squishing geometry
Ws = (5*n^2)/(4*(n-w*(2*m+b))*(w*m*(m+b)));
Gs = (Ro^2)/(Ro^3-Ri^3);

%Calculation of bending factor and bending geometry
Wb = ((2*m+b)*n)/(2*pi*w*m*(m+b));
Gb = (L^4)/(Ro*(Ro^4-Ri^4));

desired_velocity_3x1 = 0:.001:max_theoretical_velocity; %In cm/s
%Required power to weight ratio for given waves and geometry of structure and desired velocity
for i = 1:length(desired_velocity_3x1)   
    power_to_weight_3x1(i) = (Ws*Gs+Wb*Gb)*geometry_factor*desired_velocity_3x1(i);
end

figure(1)
subplot(3,1,1)
plot(desired_velocity_2x1,power_to_weight_2x1);
subplot(3,1,2)
plot(desired_velocity_2x2,power_to_weight_2x2);
subplot(3,1,3)
plot(desired_velocity_3x1,power_to_weight_3x1);

figure(2)
plot(desired_velocity_2x2,power_to_weight_2x2,desired_velocity_3x1,power_to_weight_3x1)
