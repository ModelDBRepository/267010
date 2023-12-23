clear all
close all
clc
tic 

%Parameters 
WE = 2.115; 
theta = 1;
p = 0.2;
a = 5;
tau1 = 1;
tau2 = 5; 
T = 0.5; %simulation time for each increment used to estimate the Lyapunov Exponent globally. 
xnot = 0.5:0.01:1;  %mesh, note that we used a finer xnot = -1:0.02:1, ynot = -1:0.02:1 in the manuscript, however this takes about ~4 hours to run.  
ynot = -1:0.01:1;
[mx,my] = meshgrid(xnot,ynot); 
r = WE*(mx + j*my); %eigenvalue of the weight matrix 
r = r(:);
store = zeros(length(r),1);
delay = 0.1;

%% Finite difference matrix required for method of lines. 
m = 10;
D = zeros(m,m);
%D = eye(m)*m/(2*delay);
D(m,m) = -m/delay; 
for j = 2:m-1
D(j,j-1) = m/(2*delay);
    D(j,j+1) = -m/(2*delay);    
end
D(m,m-1) = m/(delay);



%% 
parfor q = 1:length(r) %iterate of the eigenvalues 
N = 3*m;
NN = N^2 + N; 
int = rand(3*m,1);
y = zeros(NN,1); 
y(1:N,1) = int;  
y(N+1:NN) = reshape(eye(N),[N^2,1]);

tic
tot = 500;
lypc = zeros(tot,N);
cum = zeros(N,1);
%% Numerically integrate, and continuously orthogonalize, see references in manuscript for a detailed description of the algorithm. 
for nt = 1:tot

[t,ys] = ode45(@(t,y) wc_delay_lin(WE,theta,p,a,tau1,tau2,D,r(q),t,y),[0,T],y);
ys = ys(end,:); 
vec= reshape(ys(N+1:NN),[N,N]);

znorm = zeros(N,1);
% gram schmidt orthogonalization 
u = 0*vec; 
u(:,1)=vec(:,1)/sqrt(vec(:,1)'*vec(:,1));
znorm(1,1) = sqrt(vec(:,1)'*vec(:,1));
for i = 2:N 
    u(:,i) = vec(:,i); 
    for j = 1:i-1 
        u(:,i) = u(:,i) - ((u(:,j)'*u(:,i))/(u(:,j)'*u(:,j)))*u(:,j);
    end
znorm(i,1) =  sqrt(u(:,i)'*u(:,i));
u(:,i) = u(:,i)/sqrt(u(:,i)'*u(:,i));
end
cum = cum + log(znorm);

lypc(nt,:) = (cum)/(T*nt);
y(1:N) = ys(1:N); 
y(N+1:NN) = reshape(u,[N^2,1]); %use to reinitilize in the next iterate 
end
store(q) = max(lypc(end,:));
%%
toc
end

%% plot
storemat = reshape(store,[length(ynot),length(xnot)]);
imagesc(xnot,ynot,(real(storemat)),[-0.06,0]), hold on 
[surf] = contourc(xnot,ynot,real(storemat),[0,0]);
theta = 0:0.01:2*pi;
colormap('hot')
plot(cos(theta),sin(theta),'m','LineWidth',2), hold on 
xlabel('Real(r)')
ylabel('Imag(r)')
colorbar
