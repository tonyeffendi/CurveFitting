clear all
clc
close all
rng(0);

%% input
res = 0.1;
[x,y] = meshgrid(-1:res:3,0:res:6);

%% gaussian shape parameters

A = 1;
x0 = 1;
sx = 0.4;
y0 = 3;
sy = 0.9;

na = 0.32;
%f(x,y) = A exp(- (x-x0)^2/2s_x^2 + (y-y0)^2/2s_y^2)
%z = A * exp(-(0.5*(x-x0).^2/sx^2)-(0.5*(y-y0).^2/sy^2));
z = MyGaussianModel(A,x,x0,sx,y,y0,sy);
zn = z + randn(size(x))*na;

%% save to a text file
if(1)
    tmp = [x(:) y(:) zn(:)];
    txtfile = [pwd '\zn.txt'];    
    fid = fopen(txtfile,'w');
    fprintf(fid, '%f %f %f\n',tmp');    
    fclose(fid);   
end


%% plot
figure
subplot(1,2,1)
mesh(x,y,z)
xlabel('x')
ylabel('y')
title('f(x,y) = A exp(- (x-x0)^2/2s_x^2 + (y-y0)^2/2s_y^2)')
grid on

subplot(1,2,2)
mesh(x,y,zn)
xlabel('x')
ylabel('y')
title(sprintf('Noise=%0.2f\n',na))


%% non linear estimation

% param to regress: 
% A, x0, sx, y0, sy
x0 = [1 1 1 3 1];
lb = [0.1 -5 0.01 -5 0.01];
ub = [2    5  2   5  2];
options = optimoptions('lsqnonlin','Jacobian','off','TolX',1e-6,'TolFun',1e-6,...
    'Display','iter','MaxIter',30000,'MaxFunEvals',30000,'DerivativeCheck','off');
    
func = @(p) MyGaussianModel(p(1),x(:), p(2),p(3), y(:) , p(4),p(5))- zn(:) ;
%p  = lsqnonlin(func,x0,lb,ub,options)
p  = lsqnonlin(func,x0)



%% plot result
z_est = MyGaussianModel(p(1),x,p(2),p(3),y,p(4),p(5));

err = z - z_est;

figure
hold on
mesh(x,y,err)
grid on
xlabel('x')
ylabel('y')
title('err')
return;
%% g2o results
% A      = 13.6675
% x0     = 3.00104
% s_x    = 1.74378
% y0     = 9.4617
% s_y    = 2.66626

z_est2 = MyGaussianModel(A, x, x0, s_x, y,y0,s_y);

err2 = z - z_est2;

figure
hold on
mesh(x,y,err2)
grid on
xlabel('x')
ylabel('y')
title('err2 g20 diff')


%% functions
function z = MyGaussianModel(A,x,x0,sx,y,y0,sy)
  %f(x,y) = A exp(- (x-x0)^2/2s_x^2 + (y-y0)^2/2s_y^2)
  z = A * exp(-(0.5*(x-x0).^2/sx^2)-(0.5*(y-y0).^2/sy^2));
end


