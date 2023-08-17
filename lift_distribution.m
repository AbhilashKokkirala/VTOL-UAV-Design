clear;clc;close all

%% constants
rho = 1.2068;
V = 18;
S = 0.7;
CL = 0.5778;
AR = 8;
b = (AR*S)^0.5;
lambda = 0.8;
C_r = 0.3273;
C_t = C_r*lambda;

%% equivalent chord distribution
a = (S*4)/(pi*b);
bspan = linspace(0,b/2,400);
aspan = a.*((1-bspan.^2/((b/2)^2))).^0.5;
cspan = C_r-((C_r-C_t)/(b/2))*bspan;
figure(1);
plot(bspan,aspan)
hold on
plot(bspan,cspan)
hold on
plot(bspan,(cspan+aspan)/2)
c_new = (cspan+aspan)/2;
xlabel("Span")
ylabel("Chord")
legend("Elliptical","Actual","Average")
title("Schrenk's Method")

%% Shear Forces
q_dash = 0.5*rho*V^2*CL.*c_new;
L = trapz(bspan,q_dash);
p = 1.076*9.81; % point load

%% L distribution
L_dist = zeros(length(bspan),1);
for i = 2:length(bspan)
    if i<(length(bspan)*0.25)
        L_dist(i) = L-trapz(bspan(1:i),q_dash(1:i))-p;
    else
        L_dist(i) = L-trapz(bspan(1:i),q_dash(1:i));
    end
end
L_dist(1) = L-p;

%% Bending moment
m = p*bspan(length(bspan)*0.25); % due to point load 
M = trapz(bspan,q_dash.*bspan);

% M distribution
M_dist = zeros(length(bspan),1);
for i = 2:length(bspan)
    if i<(length(bspan)*0.25)
        M_dist(i) = M-m-(trapz(bspan(1:i),q_dash(1:i).*bspan(1:i)))-L_dist(i)*bspan(i);
    else
        M_dist(i) = M-trapz(bspan(1:i),q_dash(1:i).*bspan(1:i))-L_dist(i)*bspan(i);
    end
end
M_dist(1) = M-m;



%% Skin Contribution

data = importdata("NACA_4412.xlsx");

t = 0.5*10^-3;
width = zeros(35,1);
for i = 1:35
    if i<35
        width(i) = sqrt((data(i+1,1)-data(i,1))^2+(data(i+1,2)-data(i,2))^2);
    else
        width(i) = sqrt((data(1,1)-data(i,1))^2+(data(1,2)-data(i,2))^2);
    end
    
end

Area = sum(width)*t;
Area_lm = zeros(35,1);
y_c = 0;
for i = 1:35
    y_c = y_c+width(i)*t*data(i,2)/Area;
    Area_lm(i) = width(i)*t*1000;
end

figure(6);
plot(data(:,1),data(:,2),"*")
hold on
yline(y_c)

I_c = 0;
for i = 1:35
    I_c = I_c+width(i)*t*(y_c-data(i,2))^2;
end

I_skin = I_c.*cspan.^3;
figure(7)
plot(bspan,I_skin)

%% Shear Force Diagram
figure(2);
plot(bspan,L_dist,'LineWidth',0.5)
title('Shear Force Diagram')
xlim([0 bspan(end)])
xlabel("Span (m)")
ylabel("Shear Force (N)")
%saveas(gcf,'Shear_Force.png')

%% Bending Moment Diagram
figure(3);
plot(bspan,M_dist,'LineWidth',0.5)
title('Bending Moment Diagram')
xlim([0 bspan(end)])
xlabel("Span (m)")
ylabel("Bending Moment (N/m)")
%saveas(gcf,'Bending_Moment.tif')

sigma_yield = 250.3e6;
fos = 10;
sigma = sigma_yield/fos;

y_max = cspan.* 0.12.*0.5;
I = (M_dist.*(y_max'))./sigma;
t = 0.5e-3;
h = 2*(y_max- t);
b_flange = ((t.*h'.^2)./4).^(-1) .* I - h'./3;

figure(4);
plot(bspan,b_flange.*1e3,'LineWidth',1)
xlim([0 bspan(end)])
xlabel("Span (m)")
ylabel("Flange Width (mm)")
%print(gcf,'Flange_variation.png','-dpng','-r300')  

figure(5)
plot(bspan,I)
hold on
plot(bspan, I_skin)
xlabel("Span")
ylabel("Moment of Inertia")
legend("Moment of Inertia of Spar","Moment of Inertia of Skin")
title("Moment of Inertia")

