function [u, s, cx, cy, cz, rx, ry, rz] = sj_stereo(x,y,z)
% Stereo input x, y, z is coordinates in a world center caemra coordinate(wccc) of both camera's coordinate
% Output [u, s, v]is svd(Q, 0)'s output
% 
% Camera Grasshopper3(GS3-U3-50S5C-C) was used at camera model.
% Resolution : 2448x2048
% Pixels : 5.0MP
% Sensor : Sony ICX625, CCD, 2/3"
% Unit pixel size : 3.45um x 3.45um
% Focal lenth : 5.14mm
% Image plane center point : Cx = 1224, Cy= 1024 
% 
% Using pinhole camera model
%% 완성되면 지우기
% clc
% clear all
% close all

%% var, 예제를 확인하고 싶으면 인풋을 지우고 아래 주석을 없애시오.
% x = 50; 
% y = 500;
% z = 0;
xyz=[x y z 1]';
%% stereo camera characteristic
% focallength [m]
f = 0.00514;
% fixel size [m]
dv = 0.00000345;
% baseline [m]
b = 100;
% 월드 좌표계에서 본 두 카메라의 월드 좌표계 위치,  World coordinate of both cameras, [m]
% 두 카메라의 월드좌표계는 양 카메라 중심에 위치하며 카메라의 좌표계와 동일한 orientation을 가짐
% Xw = 0; Yw = 0; Zw = 0.5;
theta = 0; alpha = -pi/2; a = 0; d= 0.5;
P_Wc2Wcam = [cos(theta)*a sin(theta)*a d]';
R_Wc2Wcam = [cos(theta) -sin(theta)*cos(alpha) sin(theta)*sin(alpha);
              sin(theta) cos(theta)*cos(alpha) -cos(theta)*sin(alpha);
              0 sin(alpha) cos(alpha)];
T_Wc2Wcam = [R_Wc2Wcam P_Wc2Wcam; 0 0 0 1];
T_Wcam2Wc = [inv(R_Wc2Wcam) -inv(R_Wc2Wcam)*P_Wc2Wcam; 0 0 0 1];
% 월드 좌표계에서 본 로컬 좌표계 각 카메라의 위치 , [m]
% 카메라의 orientation은 광학축 방향이 z, 오른쪽 방향이 x, 아래쪽 방향이 y
LeftCamPos = [P_Wc2Wcam(1)-b/2,P_Wc2Wcam(2),P_Wc2Wcam(3)]';
RightCamPos = [P_Wc2Wcam(1)+b/2,P_Wc2Wcam(2),P_Wc2Wcam(3)]';
% Image plane center point
Cx = 1224; Cy= 1024;
%% Transformation Global X,Y,Z to World camera Xc,Yc,Zc (좌표 시점 변환)
% Transformation Matrix는 서로 순서가 바뀌어야 좌표 변환이 제대로 됨
% ex) World coordinate에서의 좌표를 World cam에서 알고싶을 때, Wcam2Wc를 곱해야 좌표가 알맞게 나옴.
xyz = T_Wcam2Wc*xyz;
x=xyz(1); y=xyz(2); z=xyz(3);
%% Coordinates of image plane      
% Get disparity & Jr-Cx
min_dpt = ceil(b*f/z/dv-1);
max_dpt = floor(b*f/z/dv+1);

min_Jr1 = ceil(min_dpt*(x/b-0.5));
min_Jr2 = floor(min_dpt*(x/b-0.5)+1);
max_Jr1 = ceil(max_dpt*(x/b-0.5));
max_Jr2 = floor(max_dpt*(x/b-0.5)+1);

mi1_x = [b*(0.5+min_Jr1/min_dpt) b*(0.5+(min_Jr1-1)/(min_dpt+1)) b*(0.5+(min_Jr1-1)/min_dpt) b*(0.5+min_Jr1/(min_dpt-1)) b*(0.5+min_Jr1/min_dpt)]; % pu1 pu2 pu3 pu4
mi1_z = [b*f/(min_dpt*dv) b*f/((min_dpt+1)*dv) b*f/(min_dpt*dv) b*f/((min_dpt-1)*dv) b*f/(min_dpt*dv)]; % pu1 pu2 pu3 pu4
Txz1 = T_Wc2Wcam*[mi1_x; zeros(1,length(mi1_x)); mi1_z; ones(1,length(mi1_x))];
mi1_x = Txz1(1,:);
mi1_y = Txz1(2,:);

mi2_x = [b*(0.5+min_Jr2/min_dpt) b*(0.5+(min_Jr2-1)/(min_dpt+1)) b*(0.5+(min_Jr2-1)/min_dpt) b*(0.5+min_Jr2/(min_dpt-1)) b*(0.5+min_Jr2/min_dpt)]; % pu1 pu2 pu3 pu4
mi2_z = [b*f/(min_dpt*dv) b*f/((min_dpt+1)*dv) b*f/(min_dpt*dv) b*f/((min_dpt-1)*dv) b*f/(min_dpt*dv)]; % pu1 pu2 pu3 pu4
Txz2 = T_Wc2Wcam*[mi2_x; zeros(1,length(mi2_x)); mi2_z; ones(1,length(mi2_x))];
mi2_x = Txz2(1,:);
mi2_y = Txz2(2,:);

ma1_x = [b*(0.5+max_Jr1/max_dpt) b*(0.5+(max_Jr1-1)/(max_dpt+1)) b*(0.5+(max_Jr1-1)/max_dpt) b*(0.5+max_Jr1/(max_dpt-1)) b*(0.5+max_Jr1/max_dpt)]; % pu1 pu2 pu3 pu4
ma1_z = [b*f/(max_dpt*dv) b*f/((max_dpt+1)*dv) b*f/(max_dpt*dv) b*f/((max_dpt-1)*dv) b*f/(max_dpt*dv)]; % pu1 pu2 pu3 pu4
Txz3 = T_Wc2Wcam*[ma1_x; zeros(1,length(ma1_x)); ma1_z; ones(1,length(ma1_x))];
ma1_x = Txz3(1,:);
ma1_y = Txz3(2,:);

ma2_x = [b*(0.5+max_Jr2/max_dpt) b*(0.5+(max_Jr2-1)/(max_dpt+1)) b*(0.5+(max_Jr2-1)/max_dpt) b*(0.5+max_Jr2/(max_dpt-1)) b*(0.5+max_Jr2/max_dpt)]; % pu1 pu2 pu3 pu4
ma2_z = [b*f/(max_dpt*dv) b*f/((max_dpt+1)*dv) b*f/(max_dpt*dv) b*f/((max_dpt-1)*dv) b*f/(max_dpt*dv)]; % pu1 pu2 pu3 pu4
Txz4 = T_Wc2Wcam*[ma2_x; zeros(1,length(ma2_x)); ma2_z; ones(1,length(ma2_x))];
ma2_x = Txz4(1,:);
ma2_y = Txz4(2,:);

% Wc좌표계에서 불확실성 영역 판단하기 위해 좌표 시점 변환 Wcam→Wc 
xyz = T_Wc2Wcam*xyz;
x=xyz(1); y=xyz(2); z=xyz(3);

% 판단
in = inpolygon(x,y,mi1_x,mi1_y);
if numel(x(in))==1
ret_Jrdpt=[min_Jr1 min_dpt];
xv=mi1_x;
zv=mi1_y;
end

in = inpolygon(x,y,mi2_x,mi2_y);
if numel(x(in))==1
ret_Jrdpt=[min_Jr2 min_dpt];
xv=mi2_x;
zv=mi2_y;
end

in = inpolygon(x,y,ma1_x,ma1_y);
if numel(x(in))==1
ret_Jrdpt=[max_Jr1 max_dpt];
xv=ma1_x;
zv=ma1_y;
end

in = inpolygon(x,y,ma2_x,ma2_y);
if numel(x(in))==1
ret_Jrdpt=[max_Jr2 max_dpt];
xv=ma2_x;
zv=ma2_y;
end

% Wcam좌표계에서 yl값을 얻기 위해 좌표 시점 변환 Wc→Wcam 
xyz = T_Wcam2Wc*xyz;
x=xyz(1); y=xyz(2); z=xyz(3);

% Get y
yl=floor(y*(max_dpt)/b);
min_yl = yl-1;
max_yl = yl+1;

yv1=[b*(yl)/(ret_Jrdpt(2)+1) b*(yl)/(ret_Jrdpt(2)-1) b*(yl-1)/(ret_Jrdpt(2)-1) b*(yl-1)/(ret_Jrdpt(2)+1) b*(yl)/(ret_Jrdpt(2)+1)]; % pu2, pu4, pd4, pd2, pu2
yv2=[b*(min_yl)/(ret_Jrdpt(2)+1) b*(min_yl)/(ret_Jrdpt(2)-1) b*(min_yl-1)/(ret_Jrdpt(2)) b*(min_yl-1)/(ret_Jrdpt(2)+1) b*(min_yl)/(ret_Jrdpt(2)+1)]; % pu2, pu4, pd4, pd2, pu2
yv3=[b*(max_yl)/(ret_Jrdpt(2)+1) b*(max_yl)/(ret_Jrdpt(2)-1) b*(max_yl-1)/(ret_Jrdpt(2)) b*(max_yl-1)/(ret_Jrdpt(2)+1) b*(max_yl)/(ret_Jrdpt(2)+1)]; % pu2, pu4, pd4, pd2, pu2
zv=[b*f/((ret_Jrdpt(2)+1)*dv) b*f/((ret_Jrdpt(2)-1)*dv) b*f/((ret_Jrdpt(2)-1)*dv) b*f/((ret_Jrdpt(2)+1)*dv) b*f/((ret_Jrdpt(2)+1)*dv)];

% Wc좌표를 얻기 위한 좌표 시점 변환
Tyz1 = T_Wc2Wcam*[zeros(1,length(ma2_x)); yv1; zv; ones(1,length(ma2_x))];
Tyz2 = T_Wc2Wcam*[zeros(1,length(ma2_x)); yv2; zv; ones(1,length(ma2_x))];
Tyz3 = T_Wc2Wcam*[zeros(1,length(ma2_x)); yv3; zv; ones(1,length(ma2_x))];
yv1 = Tyz1(2,:);
yv2 = Tyz2(2,:);
yv3 = Tyz3(2,:);
zv = Tyz3(3,:);

% Wc좌표계에서불확실성 영역 판단하기 위해 좌표 시점 변환 Wcam→Wc 
xyz = T_Wc2Wcam*xyz;
x=xyz(1); y=xyz(2); z=xyz(3);

in = inpolygon(y,z,yv1,zv);
if numel(y(in))==1
    yv=yv1;
    yl=yl;
end

in = inpolygon(y,z,yv2,zv);
if numel(y(in))==1
    yv=yv2;
    yl=min_yl;
end

in = inpolygon(y,z,yv3,zv);
if numel(y(in))==1
    yv=yv3;
    yl=max_yl;
end
%% 주어진 점을 포함하는 영역 확인
% figure(1)  % Get x,z
% hold on
% plot(xv,zv,'-k','LineWidth',3) % polygon
% plot(x,y,'+r','LineWidth',2)
% plot([LeftCamPos(1) RightCamPos(1)], [LeftCamPos(2) RightCamPos(2)],'ob','Linewidth',3)

% figure(1)  % Get y
% hold on
% plot(yv,zv,'LineWidth',2) % polygon
% plot(y,z,'+r','LineWidth',2)
% xlabel('distance y[m]')
% ylabel('distance z[m]')

%% Left camera & Right camera [SI:m]
Jl = ret_Jrdpt(2)+ret_Jrdpt(1);
Jr = ret_Jrdpt(1);
YL = yl;
p1 = [b*(0.5+(Jr)/(Jl-Jr)), b*YL/(Jl-Jr), b*(YL-1)/(Jl-Jr), b*f/((Jl-Jr)*dv)];
p2 = [b*(0.5+(Jr-1)/(Jl-Jr+1)), b*YL/(Jl-Jr+1), b*(YL-1)/(Jl-Jr+1), b*f/((Jl-Jr+1)*dv)];
p3 = [b*(0.5+(Jr-1)/(Jl-Jr)), b*YL/(Jl-Jr), b*(YL-1)/(Jl-Jr), b*f/((Jl-Jr)*dv)];
p4 = [b*(0.5+Jr/(Jl-Jr-1)),  b*YL/(Jl-Jr-1), b*(YL-1)/(Jl-Jr-1),  b*f/((Jl-Jr-1)*dv)];
        
Up_3D_p1= [p1(1) p1(2) p1(4) 1]';
Up_3D_p2= [p2(1) p2(2) p2(4) 1]';
Up_3D_p3= [p3(1) p3(2) p3(4) 1]';
Up_3D_p4= [p4(1) p4(2) p4(4) 1]';
Down_3D_p1= [p1(1) p1(3) p1(4) 1]';
Down_3D_p2= [p2(1) p2(3) p2(4) 1]';
Down_3D_p3= [p3(1) p3(3) p3(4) 1]';
Down_3D_p4= [p4(1) p4(3) p4(4) 1]';

% rand var <<<<<<<<
Rand_l=rand(1,1000);
Rand_t=rand(1,1000);
Rand_r=rand(1,1000);             
Rand_x = b*(0.5+(Jr+Rand_r-1)./(Jl-Jr+Rand_l-Rand_r));
Rand_z = b*f./((Jl-Jr-Rand_r+Rand_l)*dv);
Rand_y = Rand_z.*(YL-Rand_t)*dv/f;
rand_xyz = [Rand_x; Rand_y; Rand_z; ones(1,length(Rand_y))];   

% 좌표 및 rand var 시점 변환 Wcam→Wc
Wc_Up_3D_p1 = T_Wc2Wcam*Up_3D_p1;
Wc_Up_3D_p2 = T_Wc2Wcam*Up_3D_p2;
Wc_Up_3D_p3 = T_Wc2Wcam*Up_3D_p3;
Wc_Up_3D_p4 = T_Wc2Wcam*Up_3D_p4;
Wc_Down_3D_p1 = T_Wc2Wcam*Down_3D_p1;
Wc_Down_3D_p2 = T_Wc2Wcam*Down_3D_p2;
Wc_Down_3D_p3 = T_Wc2Wcam*Down_3D_p3;
Wc_Down_3D_p4 = T_Wc2Wcam*Down_3D_p4;
rand_xyz = T_Wc2Wcam*rand_xyz;
rx = rand_xyz(1,1);
ry = rand_xyz(2,1);
rz = rand_xyz(3,1);
rand_xyz = [rand_xyz(1,:)' rand_xyz(2,:)' rand_xyz(3,:)'];

%% 타원체 및 SVD 
k = convhull(rand_xyz);
xyzc=rand_xyz(k,:);
[Q, C]=MinVolEllipse(xyzc',.01);
[u s v] = svd(Q,0)
xyz_center=mean(rand_xyz);     %% xyz의 평균값 구함 -> 타원의 중심 좌표가 된다.
cx = xyz_center(1);
cy = xyz_center(2);
cz = xyz_center(3);
% 
% %% 불확실성 영역 타원
% figure(2)
% clf;
% plot_ellipse(Q,xyz_center);hold on;    %% plot_ellipse 함수를 통해 타원체를 그린다.
% plot3(rand_xyz(:,1),rand_xyz(:,2),rand_xyz(:,3),'.','markersize',20);   %% xyz 랜덤Z좌표를 점으로 표시한다.
% plot3(x,y,z,'xc','markersize',10,'linewidth',10);
% plot3(xyz_center(1),xyz_center(2),xyz_center(3),'+r','markersize',10,'linewidth',10); %% 타원 중심 좌표를 빨간색+로 표시
% title('ellipsoid modeling')
% xlabel('distance x(m)')
% ylabel('distance y(m)')
% zlabel('distance z(m)')
% set(0,'defaultfigurecolor',[1 1 1])
% hold on  %% 고정
% grid on  %% 격자 생성

%% 카메라 위치와 불확실성 영역의 위치 확인
% figure(3) 
% hold on  %% 고정
% plot3(LeftCamPos(1),LeftCamPos(2),LeftCamPos(3),'ok','markersize',10,'linewidth',8); 
% plot3(RightCamPos(1),RightCamPos(2),RightCamPos(3),'ok','markersize',10,'linewidth',8); 
% 
% plot3([Wc_Up_3D_p1(1), Wc_Up_3D_p2(1), Wc_Up_3D_p3(1), Wc_Up_3D_p4(1),Wc_Up_3D_p1(1)], ...
%       [Wc_Up_3D_p1(2), Wc_Up_3D_p2(2), Wc_Up_3D_p3(2), Wc_Up_3D_p4(2),Wc_Up_3D_p1(2)], ...
%       [Wc_Up_3D_p1(3), Wc_Up_3D_p2(3), Wc_Up_3D_p3(3), Wc_Up_3D_p4(3),Wc_Up_3D_p1(3)], '-r','linewidth',1)
% plot3([Wc_Down_3D_p1(1), Wc_Down_3D_p2(1), Wc_Down_3D_p3(1), Wc_Down_3D_p4(1),Wc_Down_3D_p1(1)], ...
%       [Wc_Down_3D_p1(2), Wc_Down_3D_p2(2), Wc_Down_3D_p3(2), Wc_Down_3D_p4(2),Wc_Down_3D_p1(2)], ...
%       [Wc_Down_3D_p1(3), Wc_Down_3D_p2(3), Wc_Down_3D_p3(3), Wc_Down_3D_p4(3),Wc_Down_3D_p1(3)], '-r','linewidth',1)
% 
% plot3([Wc_Up_3D_p1(1), Wc_Down_3D_p1(1)], [Wc_Up_3D_p1(2), Wc_Down_3D_p1(2)], [Wc_Up_3D_p1(3), Wc_Down_3D_p1(3)], '-k','linewidth',1)
% plot3([Wc_Up_3D_p2(1), Wc_Down_3D_p2(1)], [Wc_Up_3D_p2(2), Wc_Down_3D_p2(2)], [Wc_Up_3D_p2(3), Wc_Down_3D_p2(3)], '-k','linewidth',1)
% plot3([Wc_Up_3D_p3(1), Wc_Down_3D_p3(1)], [Wc_Up_3D_p3(2), Wc_Down_3D_p3(2)], [Wc_Up_3D_p3(3), Wc_Down_3D_p3(3)], '-k','linewidth',1)
% plot3([Wc_Up_3D_p4(1), Wc_Down_3D_p4(1)], [Wc_Up_3D_p4(2), Wc_Down_3D_p4(2)], [Wc_Up_3D_p4(3), Wc_Down_3D_p4(3)], '-k','linewidth',1)
% plot_ellipse(Q,xyz_center);hold on;    %% plot_ellipse 함수를 통해 타원체를 그린다.
% plot3(rand_xyz(:,1),rand_xyz(:,2),rand_xyz(:,3),'.','markersize',5);   %% xyz 랜덤Z좌표를 점으로 표시한다.
% plot3(xyz_center(1),xyz_center(2),xyz_center(3),'+c','markersize',20); %% 타원 중심 좌표를 빨간색+로 표시
% title('ellipsoid modeling')
% xlabel('distance x(m)')
% ylabel('distance z(m)')
% zlabel('distance y(m)')
end


function [A , c] = MinVolEllipse(P, tolerance)
% [A , c] = MinVolEllipse(P, tolerance)
% Finds the minimum volume enclsing ellipsoid (MVEE) of a set of data
% points stored in matrix P. The following optimization problem is solved: 
%
% minimize       log(det(A))
% subject to     (P_i - c)' * A * (P_i - c) <= 1
%                
% in variables A and c, where P_i is the i-th column of the matrix P. 
% The solver is based on Khachiyan Algorithm, and the final solution 
% is different from the optimal value by the pre-spesified amount of 'tolerance'.
%
% inputs:
%---------
% P : (d x N) dimnesional matrix containing N points in R^d.
% tolerance : error in the solution with respect to the optimal value.
%
% outputs:
%---------
% A : (d x d) matrix of the ellipse equation in the 'center form': 
% (x-c)' * A * (x-c) = 1 
% c : 'd' dimensional vector as the center of the ellipse. 
% 
% example:
% --------
%      P = rand(5,100);
%      [A, c] = MinVolEllipse(P, .01)
%
%      To reduce the computation time, work with the boundary points only:
%      
%      K = convhulln(P');  
%      K = unique(K(:));  
%      Q = P(:,K);
%      [A, c] = MinVolEllipse(Q, .01)
%
%
% Nima Moshtagh (nima@seas.upenn.edu)
% University of Pennsylvania
%
% December 2005
% UPDATE: Jan 2009

%%%%%%%%%%%%%%%%%%%%% Solving the Dual problem%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% ---------------------------------
% data points 
% -----------------------------------
[d N] = size(P);

Q = zeros(d+1,N);
Q(1:d,:) = P(1:d,1:N);
Q(d+1,:) = ones(1,N);

% initializations
% -----------------------------------
count = 1;
err = 1;
u = (1/N) * ones(N,1);          % 1st iteration

% Khachiyan Algorithm
% -----------------------------------
while err > tolerance
    X = Q * diag(u) * Q';       % X = \sum_i ( u_i * q_i * q_i')  is a (d+1)x(d+1) matrix
    M = diag(Q' * inv(X) * Q);  % M the diagonal vector of an NxN matrix
    [maximum j] = max(M);
    step_size = (maximum - d -1)/((d+1)*(maximum-1));
    new_u = (1 - step_size)*u ;
    new_u(j) = new_u(j) + step_size;
    count = count + 1;
    err = norm(new_u - u);
    u = new_u;
end

%%%%%%%%%%%%%%%%%%% Computing the Ellipse parameters%%%%%%%%%%%%%%%%%%%%%%
% Finds the ellipse equation in the 'center form': 
% (x-c)' * A * (x-c) = 1
% It computes a dxd matrix 'A' and a d dimensional vector 'c' as the center
% of the ellipse. 

U = diag(u);

% the A matrix for the ellipse
% --------------------------------------------
A = (1/d) * inv(P * U * P' - (P * u)*(P*u)' );
A = inv(A);

% center of the ellipse 
% --------------------------------------------
c = P * u;
end
