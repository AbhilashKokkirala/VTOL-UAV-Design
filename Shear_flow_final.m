clear;clc;close all;
%% To extend the no of points
load("Lift_dist.mat");
data = importdata("NACA_4412.xlsx");
ind = find(data(:,1)==0);
upper_points = data(1:ind,:);
lower_points = data(ind+1:end,:);
x_upper = 1:-0.05:0.05;
x_lower = 0.05:0.05:1;
x_upper_new = zeros(length(x_upper)+1,1);
x_upper_new(1:length(x_upper),1) = x_upper;
x_upper_new(end) = 0;
x_upper_new(x_upper_new == 0.3) = [];
x_lower(x_lower == 0.3) = [];
upper_points_inter = interp1(upper_points(:,1),upper_points(:,2),x_upper_new,'spline');
lower_points_inter = interp1(lower_points(:,1),lower_points(:,2),x_lower','spline');
data_inter = zeros(length(upper_points_inter)+length(lower_points_inter),2);
data_inter(1:length(upper_points_inter),1) = x_upper_new;
data_inter(1:length(upper_points_inter),2) = upper_points_inter;
data_inter(length(upper_points_inter)+1:length(upper_points_inter) + length(lower_points_inter),1) = x_lower';
data_inter(length(upper_points_inter)+1:length(upper_points_inter) + length(lower_points_inter),2) = lower_points_inter;

data = data_inter;

% x = [1,0.6,0.075,0,0.075,0.6];
% y = [-0.0013, 0.0814, 0.0576, 0, -0.0274, -0.01];
% data = zeros(length(x),2);
% data(:,1) = x;
% data(:,2) = y;

%% Shear flow without spar

t = 0.5*1e-3;
V_y = L_dist;
c = cspan;
ix = length(data);
iy = length(cspan);
y = c.*data(:,2);
x = c.*data(:,1);
x_c = zeros(1,iy);
y_c = zeros(1,iy);
A_bar = zeros(1,iy);
width = zeros(ix,1);
for i = 1:ix % Width,(x,y) of centriod, A_bar
    if i<ix
        width(i) = sqrt((data(i+1,1)-data(i,1))^2+(data(i+1,2)-data(i,2))^2);
    else
        width(i) = sqrt((data(1,1)-data(i,1))^2+(data(1,2)-data(i,2))^2);
    end
    
end
width = c.*width;
Area = zeros(1,iy);
for i = 1:iy
    Area(i) = sum(width(:,i))*t;
end
for i =1:ix
        x_c = x_c+width(i,:).*t.*x(i,:)./(Area);
        y_c = y_c+width(i,:).*t.*y(i,:)./(Area);
        A_bar = A_bar+abs(width(i,:).*y(i,:));
end

I_c_xx = zeros(1,iy);
I_c_yy = zeros(1,iy);
I_c_xy = zeros(1,iy);

I_x = width.*t^3./12;
I_y = (width).^3.*t./12;
I_u = zeros(ix,iy);
I_v = zeros(ix,iy);
I_uv = zeros(ix,iy);
for i = 1:ix
    if i<ix
        I_u(i,:) = (I_x(i)+I_y(i))./2+(I_x(i)-I_y(i))./2.*((x(i+1,:)-x(i,:)).^2-(y(i+1,:)-y(i,:)).^2)./((x(i+1,:)-x(i,:)).^2+(y(i+1,:)-y(i,:)).^2);
        I_c_xx = I_c_xx+I_u(i)+width(i,:)*t.*(y_c-(y(i,:)+y(i+1,:))./2).^2;
        
        I_v(i,:) = (I_x(i)+I_y(i))./2-(I_x(i)-I_y(i))./2.*((x(i+1,:)-x(i,:)).^2-(y(i+1,:)-y(i,:)).^2)./((x(i+1,:)-x(i,:)).^2+(y(i+1,:)-y(i,:)).^2);
        I_c_yy = I_c_yy+I_v(i)+width(i,:)*t.*(x_c-(x(i,:)+x(i+1,:))./2).^2;
        
        I_uv(i,:) = (I_x(i)-I_y(i)).*((x(i+1,:)-x(i,:)).*(y(i+1,:)-y(i,:)))./((x(i+1,:)-x(i,:)).^2+(y(i+1,:)-y(i,:)).^2);
        I_c_xy = I_c_xy+I_uv(i)+width(i,:)*t.*(x_c-(x(i,:)+x(i+1,:))./2).*(y_c-(y(i,:)+y(i+1,:))./2);
        
    else
        I_u(i,:) = (I_x(i)+I_y(i))./2+(I_x(i)-I_y(i))./2.*((x(1,:)-x(i,:)).^2-(y(1,:)-y(i,:)).^2)./((x(1,:)-x(i,:)).^2+(y(1,:)-y(i,:)).^2);
        I_c_xx = I_c_xx+I_u(i)+width(i,:)*t.*(y_c-(y(1,:)+y(i,:))./2).^2;
        
        I_v(i,:) = (I_x(i)+I_y(i))./2-(I_x(i)-I_y(i))./2.*((x(1,:)-x(i,:)).^2-(y(1,:)-y(i,:)).^2)./((x(1,:)-x(i,:)).^2+(y(1,:)-y(i,:)).^2);
        I_c_yy = I_c_yy+I_v(i)+width(i,:)*t.*(x_c-(x(1,:)+x(i,:))./2).^2;
        
        I_uv(i,:) = (I_x(i)-I_y(i)).*((x(1,:)-x(i,:)).*(y(1,:)-y(i,:)))./((x(1,:)-x(i,:)).^2+(y(1,:)-y(i,:)).^2);
        I_c_xy = I_c_xy+I_uv(i)+width(i,:)*t.*(x_c-(x(1,:)+x(i,:))./2).*(y_c-(y(i,:)+y(1,:))./2);
    end
end

I_skin_xx = I_c_xx;
I_skin_yy = I_c_yy;
I_skin_xy = I_c_xy;

I_skin_eq = (I_skin_xx.*I_skin_yy-I_skin_xy.^2);
K_x = I_skin_xx./I_skin_eq;
K_y = I_skin_yy./I_skin_eq;
K_xy = I_skin_xy./I_skin_eq;

sigma_skin =  ( (y-y_c) .* ( cspan.*M_dist'.*K_y) )-( (x-x_c) .* (cspan.*M_dist'.*K_xy) );

area_lm = zeros(ix,iy);

for i = 1:ix % limped mass area
    if i == 1
        area_lm(i,:) = t.*width(ix,:)/6.*(2+sigma_skin(ix,:)./sigma_skin(i,:))+t*width(i,:)/6.*(2+sigma_skin(i+1,:)./sigma_skin(i,:));
    elseif i == ix
        area_lm(i,:) = t.*width(i-1,:)/6.*(2+sigma_skin(i-1,:)./sigma_skin(i,:))+t*width(i,:)/6.*(2+sigma_skin(1,:)./sigma_skin(i,:));
    else
        area_lm(i,:) = t.*width(i-1,:)/6.*(2+sigma_skin(i-1,:)./sigma_skin(i,:))+t*width(i,:)/6.*(2+sigma_skin(i+1,:)./sigma_skin(i,:));
    end
end

%q_dash and M
Q_x = zeros(ix,iy);
Q_y = zeros(ix,iy);
q_dash = zeros(ix,iy);
M = zeros(1,iy);

for i = 1:ix-1
    if i >1
        Q_x(i,:) = Q_x(i-1,:)+area_lm(i+1,:).*(y(i+1,:)-y_c);
        Q_y(i,:) = Q_y(i-1,:)+area_lm(i+1,:).*(x(i+1,:)-x_c);
    else
        Q_x(1,:) = area_lm(2,:).*(y(2,:)-y_c);
        Q_y(1,:) = area_lm(2,:).*(x(2,:)-x_c);
    end
    q_dash(i+1,:) = ( -( Q_x(i,:).*( V_y'.*K_y) )+( Q_y(i,:).*(V_y'.*K_xy) ) );%( -( Q_x(i,:).*( V_y')./I_c_xx ));%
    if i < ix-1
        M = M-q_dash(i+1,:).*width(i+1,:).*abs((y(i+2,:)-y(i+1,:)).*x_c+(x(i+1,:)-x(i+2,:)).*y_c+(x(i+2,:)-x(i+1,:)).*y(i+1,:)-(y(i+2,:)-y(i+1,:)).*x(i+1,:))./(sqrt((y(i+2,:)-y(i+1,:)).^2+(x(i+1,:)-x(i+2,:)).^2));
    else
        M = M-q_dash(i+1,:).*width(i+1,:).*abs((y(1,:)-y(i+1,:)).*x_c+(x(i+1,:)-x(1,:)).*y_c+(x(1,:)-x(i+1,:)).*y(i+1,:)-(y(1,:)-y(i+1,:)).*x(i+1,:))./(sqrt((y(1,:)-y(i+1,:)).^2+(x(i+1,:)-x(1,:)).^2));
    end
end

q_o = M./(2.*A_bar); 
q = q_o.*ones(ix,1)+q_dash;

figure(1)
scatter3(reshape(bspan.*ones(ix,1),[],1),reshape(x,[],1),reshape(y,[],1),[],reshape(q,[],1),'fill')
colorbar

figure(2)
j=1;
fprintf("%d",V_y(j))
scatter(x(:,j),y(:,j),[],q(:,j),'fill')
ylim([-0.05,0.05])
title ("Shear flow")
colorbar

Shear_flow_skin = zeros(ix,iy,3);
Shear_flow_skin(:,:,1) = x;
Shear_flow_skin(:,:,2) = y;
Shear_flow_skin(:,:,3) = q;

%% Shear flow with spar
%% Calculation of I
a = t;
H = 0.12*cspan;
h = H - 2*t;
b_flange = 40e-3; 

Area_spar = a.*h + (2*b_flange).*(H-h);
x_spar = 0.3.*cspan;
x_c_new = (x_c.*Area + x_spar.*Area_spar)./(Area+Area_spar);

I_spar_xx = a.*(h.^3)./12 + (b_flange/12).*(H.^3-h.^3); %%I_spar_xx = a.*(h.*3)./12 + (b_flange/12).*(H.^3-h.^3);
I_spar_yy = (a^3).*h./12 + ((b_flange^3)/12).*(H-h) + ((x_c_new-x_spar).^2.*Area_spar).*cspan.^2;

I_skin_yy_new = I_skin_yy + ((x_c_new-x_c).^2.*Area_spar).*cspan.^2;

I_total_xx = I_skin_xx + I_spar_xx;
I_total_yy = I_skin_yy_new + I_spar_yy;
I_total_xy = I_skin_xy;
I_total_eq = (I_total_xx.*I_total_yy-I_total_xy.^2);

K_x = I_total_xx./I_total_eq;
K_y = I_total_yy./I_total_eq;
K_xy = I_total_xy./I_total_eq;

%% Spliting the loop
%spar_points = [[0.3,0.1];[0.3,0.08];[0.3,0.06];[0.3,0.04];[0.3,0.02];[0.3,0];[0.3,-0.02]];
spar_points = [[0.3,0.0976];[0.3,0.04];[0.3,-0.0226]];
j=0;
k=0;
change = "False";
for i = 1:ix % Spliting into r and l loop

    if data(i,1) >= 0.3
        if change == "False"
            k = k+1;
            right_loop_2(k,1) = data(i,1);
            right_loop_2(k,2) = data(i,2);
        else
            k = k+1;
            right_loop_1(k,1) = data(i,1);
            right_loop_1(k,2) = data(i,2);
        end
    end
    
    if data(i,1) <= 0.3
        j=j+1;
        change = "True";
        k=0;
        left_loop(j,1) = data(i,1);
        left_loop(j,2) = data(i,2);
    end
    
end


for i = 1:length(spar_points)
    left_loop(i+j,1) = spar_points(length(spar_points)+1-i,1);
    left_loop(i+j,2) = spar_points(length(spar_points)+1-i,2);
end

right_loop = zeros(length(right_loop_1)+length(right_loop_2)+2,2);
right_loop(1,:) = spar_points(end,:);
right_loop(2:1+length(right_loop_1),:) = right_loop_1;
%right_loop(2,:) = right_loop_1;
right_loop(length(right_loop_1)+2:length(right_loop_1)+length(right_loop_2)+1,:) = right_loop_2;
right_loop(end,:) = spar_points(1,:);

figure(11);
plot(left_loop(:,1),left_loop(:,2),"*")
hold on
plot(spar_points(:,1),spar_points(:,2),"o")

figure(12);
plot(right_loop(:,1),right_loop(:,2),"*")
hold on
plot(spar_points(:,1),spar_points(:,2),"o")

%%
ix_l = length(left_loop);
width_lloop = zeros(ix_l,1); %length of each segment
for i = 1:ix_l
    if i<ix_l
        width_lloop(i) = sqrt((left_loop(i+1,1)-left_loop(i,1))^2+(left_loop(i+1,2)-left_loop(i,2))^2);
    else
        width_lloop(i) = sqrt((left_loop(1,1)-left_loop(i,1))^2+(left_loop(1,2)-left_loop(i,2))^2);
    end
    
end
Area_lloop = sum(width_lloop)*t; % Material area of left loop
y = c.*left_loop(:,2);
x = c.*left_loop(:,1);

x_c_lloop = zeros(1,iy);
y_c_lloop = zeros(1,iy);
A_1bar = zeros(1,iy); % Area enclosed by left loop
for i = 1:ix_l
    if x(i)>c*0.3-0.02
        x_c_lloop = x_c_lloop+c.*width_lloop(i)*2*t.*x(i,:)./(c.*Area_lloop);
        y_c_lloop = y_c_lloop+c.*width_lloop(i)*2*t.*y(i,:)./(c.*Area_lloop);
    else
        x_c_lloop = x_c_lloop+c.*width_lloop(i)*t.*x(i,:)./(c.*Area_lloop);
        y_c_lloop = y_c_lloop+c.*width_lloop(i)*t.*y(i,:)./(c.*Area_lloop);
    end
    if (i <ix_l-length(spar_points)) || (i==ix_l)
        A_1bar = A_1bar+abs(c.*width_lloop(i).*y(i,:));
    end
end


ix_r = length(right_loop);
width_rloop = zeros(ix_r,1);
for i = 1:ix_r-1 % Last lumped area has no area
    width_rloop(i) = sqrt((right_loop(i+1,1)-right_loop(i,1))^2+(right_loop(i+1,2)-right_loop(i,2))^2);
end
Area_rloop = sum(width_rloop)*t; %Material area of right loop

y = c.*right_loop(:,2);
x = c.*right_loop(:,1);

x_c_rloop = zeros(1,iy);
y_c_rloop = zeros(1,iy);
A_2bar = zeros(1,iy); % Area enclosed by right loop
for i = 1:ix_r
    if x(i)<c*0.3+0.02
        x_c_rloop = x_c_rloop+c.*width_rloop(i)*2*t.*x(i,:)./(c.*Area_rloop);
        y_c_rloop = y_c_rloop+c.*width_rloop(i)*2*t.*y(i,:)./(c.*Area_rloop);
    else
        x_c_rloop = x_c_rloop+c.*width_rloop(i)*t.*x(i,:)./(c.*Area_rloop);
        y_c_rloop = y_c_rloop+c.*width_rloop(i)*t.*y(i,:)./(c.*Area_rloop);
    end
    A_2bar = A_2bar+abs(c.*width_rloop(i).*y(i,:));
end

x_c = (x_c_lloop.*Area_lloop + x_c_rloop.*Area_rloop)./(Area_lloop+Area_rloop);
y_c = (y_c_lloop.*Area_lloop + y_c_rloop.*Area_rloop)./(Area_lloop+Area_rloop);
%%

x_leftloop = left_loop(:,1);
y_leftloop = left_loop(:,2);
x_rightloop = right_loop(:,1);
y_rightloop = right_loop(:,2);

sigma_total_leftloop =  ( (y_leftloop-y_c).*( cspan.*M_dist'.*K_y ) )-( (x_leftloop-x_c).*(cspan.*M_dist'.*K_xy) ); % Verified
sigma_total_rightloop =  ( (y_rightloop-y_c).*( cspan.*M_dist'.*K_y ) )-( (x_rightloop-x_c).*(cspan.*M_dist'.*K_xy) );

save("sigma_total_leftloop")
save("sigma_total_rightloop")

%%
index_left = find ((left_loop(:,1)>0.3-1e-6));
top_left = length(left_loop);
bottom_left = index_left(1);

index_right = find ((right_loop(:,1)<0.3+1e-6));
top_right = index_right(2)-1;
bottom_right = 2;

width_top = sqrt((left_loop(top_left,1)-right_loop(top_right,1))^2+(left_loop(top_left,2)-right_loop(top_right,2))^2);
width_bottom = sqrt((left_loop(bottom_left,1)-right_loop(bottom_right,1))^2+(left_loop(bottom_left,2)-right_loop(bottom_right,2))^2);

y = c.*left_loop(:,2);
x = c.*left_loop(:,1);

area_lm_lloop = zeros(ix_l,iy);
for i = 1:ix_l % Lumping area for left loop
    if i == top_left
        area_lm_lloop(i,:) = t.*c.*width_lloop(i-1)/6.*(2+sigma_total_leftloop(i-1,:)./sigma_total_leftloop(i,:))+t.*c.*width_lloop(i)/6.*(2+sigma_total_leftloop(1,:)./sigma_total_leftloop(i,:))+t*c.*width_top/6.*(2+sigma_total_rightloop(top_right,:)./sigma_total_leftloop(i,:));
    elseif i == bottom_left
        area_lm_lloop(i,:) = t.*c.*width_lloop(i-1)/6.*(2+sigma_total_leftloop(i-1,:)./sigma_total_leftloop(i,:))+t.*c.*width_lloop(i)/6.*(2+sigma_total_leftloop(i+1,:)./sigma_total_leftloop(i,:))+t*c.*width_bottom/6.*(2+sigma_total_rightloop(bottom_right,:)./sigma_total_leftloop(i,:));
    elseif i == 1
        area_lm_lloop(i,:) = t.*c.*width_lloop(ix_l)/6.*(2+sigma_total_leftloop(ix_l,:)./sigma_total_leftloop(i,:))+t*c.*width_lloop(i)/6.*(2+sigma_total_leftloop(i+1,:)./sigma_total_leftloop(i,:));
    else
        area_lm_lloop(i,:) = t.*c.*width_lloop(i-1)/6.*(2+sigma_total_leftloop(i-1,:)./sigma_total_leftloop(i,:))+t*c.*width_lloop(i)/6.*(2+sigma_total_leftloop(i+1,:)./sigma_total_leftloop(i,:));
    end
end

Q_x = zeros(ix_l,iy);
Q_y = zeros(ix_l,iy);
q_dash = zeros(ix_l,iy);

for i = 1:ix_l-1 % Cut at the last skin portion
    if i >1
        Q_x(i,:) = Q_x(i-1,:)+area_lm_lloop(i+1,:).*(y(i+1,:)-y_c_lloop);
        Q_y(i,:) = Q_y(i-1,:)+area_lm_lloop(i+1,:).*(x(i+1,:)-x_c_lloop);
    else
        Q_x(1,:) = area_lm_lloop(1,:).*(y(1,:)-y_c_lloop);
        Q_y(1,:) = area_lm_lloop(1,:).*(x(1,:)-x_c_lloop);
    end
    q_dash(i,:) = ( -( Q_x(i,:).*( V_y'.*K_y) )+( Q_y(i,:).*(V_y'.*K_xy) ) );%( -( Q_x(i,:).*( V_y')./I_c_xx ));% Verified
end

figure(3)
scatter3(reshape(bspan.*ones(ix_l,1),[],1),reshape(x,[],1),reshape(y,[],1),[],reshape(q_dash,[],1),'fill')
colorbar

figure(4)
j=100;
scatter(x(:,j),y(:,j),[],q_dash(:,j),'fill')
ylim([-0.05,0.05])
title ("Shear flow")
colorbar

Shear_flow_lloop = zeros(ix_l,iy,3);
Shear_flow_lloop(:,:,1) = x;
Shear_flow_lloop(:,:,2) = y;
Shear_flow_lloop(:,:,3) = q_dash;
%%

y = c.*right_loop(:,2);
x = c.*right_loop(:,1);

area_lm_rloop = zeros(ix_r,iy);
for i = 1:ix_r
    if i == ix_r
        area_lm_rloop(i,:) = area_lm_lloop(top_left,:);
    elseif i == 1
        area_lm_rloop(i,:) = area_lm_lloop(bottom_left,:);
    else
        area_lm_rloop(i,:) = t.*c.*width_rloop(i-1)/6.*(2+sigma_total_rightloop(i-1,:)./sigma_total_rightloop(i,:))+t*c.*width_rloop(i)/6.*(2+sigma_total_rightloop(i+1,:)./sigma_total_rightloop(i,:));
    end
end

Q_x = zeros(ix_r,iy);
Q_y = zeros(ix_r,iy);
q_dash = zeros(ix_r,iy);

for i = 1:ix_r-1 % cut at the first skin portion
    if i >1
        Q_x(i,:) = Q_x(i-1,:)+area_lm_rloop(i+1,:).*(y(i+1,:)-y_c_rloop);
        Q_y(i,:) = Q_y(i-1,:)+area_lm_rloop(i+1,:).*(x(i+1,:)-x_c_rloop);
    else
        Q_x(1,:) = area_lm_rloop(2,:).*(y(2,:)-y_c_rloop);
        Q_y(1,:) = area_lm_rloop(2,:).*(x(2,:)-x_c_rloop);
    end
    q_dash(i+1,:) = ( -( Q_x(i,:).*( V_y'.*K_y) )+( Q_y(i,:).*(V_y'.*K_xy) ) );%( -( Q_x(i,:).*( V_y')./I_c_xx ));%
end

figure(5)
scatter3(reshape(bspan.*ones(ix_r,1),[],1),reshape(x,[],1),reshape(y,[],1),[],reshape(q_dash,[],1),'fill')
colorbar

figure(6)
j=100;
scatter(x(:,j),y(:,j),[],q_dash(:,j),'fill')
ylim([-0.05,0.05])
title ("Shear flow")
colorbar

Shear_flow_rloop = zeros(ix_r,iy,3);
Shear_flow_rloop(:,:,1) = x;
Shear_flow_rloop(:,:,2) = y;
Shear_flow_rloop(:,:,3) = q_dash;
%% Finding q1 and q2

%left loop
M = zeros(1,iy-1); % Moment about centroid
q_dash = Shear_flow_lloop(:,1:iy-1,3);
y = Shear_flow_lloop(:,1:iy-1,2);
x = Shear_flow_lloop(:,1:iy-1,1);
for i = 1:ix_l
    if i<ix_l
        M = M-q_dash(i,:).*c(1:iy-1).*width_lloop(i).*abs((y(i+1,:)-y(i,:)).*x_c(1:iy-1)+(x(i,:)-x(i+1,:)).*y_c(1:iy-1)+(x(i+1,:)-x(i,:)).*y(i,:)-(y(i+1,:)-y(i,:)).*x(i,:))./(sqrt((y(i+1,:)-y(i,:)).^2+(x(i,:)-x(i+1,:)).^2));
    else
        M = M-q_dash(i,:).*c(1:iy-1).*width_lloop(i).*abs((y(1,:)-y(i,:)).*x_c(1:iy-1)+(x(i,:)-x(1,:)).*y_c(1:iy-1)+(x(1,:)-x(i,:)).*y(i,:)-(y(1,:)-y(i,:)).*x(i,:))./(sqrt((y(1,:)-y(i,:)).^2+(x(i,:)-x(1,:)).^2));
    end
end

%right loop
q_dash = Shear_flow_rloop(:,1:iy-1,3);
y = Shear_flow_rloop(:,1:iy-1,2);
x = Shear_flow_rloop(:,1:iy-1,1);
for i = 1:ix_r-1 % anyway last width is zero
    M = M-q_dash(i,:).*c(1:iy-1).*width_rloop(i).*abs((y(i+1,:)-y(i,:)).*x_c(1:iy-1)+(x(i,:)-x(i+1,:)).*y_c(1:iy-1)+(x(i+1,:)-x(i,:)).*y(i,:)-(y(i+1,:)-y(i,:)).*x(i,:))./(sqrt((y(i+1,:)-y(i,:)).^2+(x(i,:)-x(i+1,:)).^2));
end

q_1 = sym("q_1",[1 iy-1]); % q_1 - integration constant for left loop
q_2 = sym("q_2",[1 iy-1]); % q_2 - integration constant for right loop

theta_lloop = zeros(1,iy-1); % dtheta/dx for left loop
theta_rloop = zeros(1,iy-1); % dtheta/dx for right loop

for i = 1:ix_l-length(spar_points)
    if left_loop(i,1)<c*0.3-0.02
        theta_lloop = theta_lloop + (c(1:iy-1).*width_lloop(i).*(q_1 + Shear_flow_lloop(i,1:iy-1,3)))./(t.*A_1bar(1:iy-1));
    else
        theta_lloop = theta_lloop + (c(1:iy-1).*width_lloop(i).*(q_1 + Shear_flow_lloop(i,1:iy-1,3)))./(2*t.*A_1bar(1:iy-1));
    end
end
theta_lloop = theta_lloop + (c(1:iy-1).*width_lloop(end).*(q_1 + Shear_flow_lloop(end,1:iy-1,3)))./(t.*A_1bar(1:iy-1)); % Final skin is not in spar

for i = 1:ix_r
    if right_loop(i,1)>c*0.3+0.02
        theta_rloop = theta_rloop + c(1:iy-1).*width_rloop(i).*(q_2 + Shear_flow_rloop(i,1:iy-1,3))./(t.*A_2bar(1:iy-1));
    else
        theta_rloop = theta_rloop + c(1:iy-1).*width_rloop(i).*(q_2 + Shear_flow_rloop(i,1:iy-1,3))./(2*t.*A_2bar(1:iy-1));
    end
end

for i = ix_l-length(spar_points)+1:ix_l-1
    theta_lloop = theta_lloop + (c(1:iy-1).*width_lloop(i).*(q_1 - q_2 + Shear_flow_lloop(i,1:iy-1,3)))./(t.*A_1bar(1:iy-1));
    theta_rloop = theta_rloop + c(1:iy-1).*width_rloop(i).*(q_2 - q_1 + Shear_flow_rloop(i,1:iy-1,3))./(t.*A_2bar(1:iy-1));
end

eqn = [2.*A_1bar(1:iy-1).*q_1+2.*A_2bar(1:iy-1).*q_2 == M, theta_lloop == theta_rloop]; % Solving net moment = 0 and theta of both loops are equal

S = solve(eqn,[q_1 q_2]);

q_dash = Shear_flow_lloop(:,1:iy-1,3);
y_lloop = Shear_flow_lloop(:,1:iy-1,2);
x_lloop = Shear_flow_lloop(:,1:iy-1,1);
q_lloop = zeros(ix_l,iy-1);
for i = 1:iy-1
    q_lloop(:,i) = q_dash(:,i) + eval(strcat("S.q_1",string(i))).*ones(ix_l,1);
end

q_dash = Shear_flow_rloop(1:ix_r,1:iy-1,3);
y_rloop = Shear_flow_rloop(1:ix_r,1:iy-1,2);
x_rloop = Shear_flow_rloop(1:ix_r,1:iy-1,1);
q_rloop = zeros(ix_r,iy-1);
for i = 1:iy-1
    q_rloop(:,i) = q_dash(:,i) + eval(strcat("S.q_2",string(i)));
end
load("shearflow_from_torsion.mat")

q_spar = zeros(ix_l+ix_r,iy-1,3);
width_loop = zeros(ix_l+ix_r,1);
q_spar(1:ix_l,:,1) = x_lloop;
q_spar(1:ix_l,:,2) = y_lloop;
q_spar(1:ix_l,:,3) = q_lloop ;%- ones(ix_l,1)*q_left_torsion(1:iy-1);
q_spar(ix_l+1:end,:,1) = x_rloop;
q_spar(ix_l+1:end,:,2) = y_rloop;
q_spar(ix_l+1:end,:,3) = q_rloop ;%- ones(ix_r,1)*q_right_torsion(1:iy-1);

width_loop(1:ix_l,1) = width_lloop;
width_loop(1:ix_r,1) = width_rloop;

x = q_spar(:,:,1);
y = q_spar(:,:,2);
q = q_spar(:,:,3);

figure(7)
scatter3(reshape(bspan(1:iy-1).*ones(ix_l+ix_r,1),[],1),reshape(x,[],1),reshape(y,[],1),[],reshape(q,[],1),'fill')
colorbar

figure(8)
j=100;
scatter(x(:,j),y(:,j),[],q(:,j),'fill')
ylim([-0.05,0.05])
title ("Shear flow")
colorbar

figure(9)
j=100;
scatter3(x(:,j),y(:,j),q(:,j))
title ("Shear flow")

%%
M = zeros(1,iy-1);
x = q_spar(:,:,1);
y = q_spar(:,:,2);
q = q_spar(:,:,3);
for i = 1:ix_l+ix_r-1 % anyway last width is zero
    M = M-q(i,:).*c(1:iy-1).*width_loop(i).*abs((y(i+1,:)-y(i,:)).*x_c(1:iy-1)+(x(i,:)-x(i+1,:)).*y_c(1:iy-1)+(x(i+1,:)-x(i,:)).*y(i,:)-(y(i+1,:)-y(i,:)).*x(i,:))./(sqrt((y(i+1,:)-y(i,:)).^2+(x(i,:)-x(i+1,:)).^2));
end

e_x = M./(V_y(1:iy-1)'.*c(1:iy-1));

CP = 0.327.*c(1:iy-1);
e_x = x_c(1:iy-1)+e_x;
e = e_x - CP;

T_lift = V_y(1:iy-1)'.*e;

q_o_lift = T_lift./(2.*A_bar(1:iy-1));

T_boom = 1.076*9.81*0.509647;% Torque = m*g*distance

q_o_boom = T_boom./(2.*A_bar(1:iy-1));
%%
per_left = c.*sum(width_lloop(1:ix_l-length(spar_points)));
per_left = per_left+c.*width_lloop(end);
per_spar = c.*sum(width_lloop(ix_l-length(spar_points)+1:end-1));
per_right = c.*sum(width_rloop(1:end));


