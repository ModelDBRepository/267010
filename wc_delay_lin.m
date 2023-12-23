function dy=wc_delay_lin(WE,theta,p,a,tau1,tau2,D,r,t,y)

m = size(D,1);
%% should be 3m + (3m)^2 total equations here! 
e =y(1:m); 
i=y(m+1:2*m);
wi=y(2*m+1:3*m);

Y = reshape(y(3*m+1:3*m+(3*m)^2),[3*m,3*m]);
dy = zeros(3*m+(3*m)^2,1);
%% WC equation for the single recurrently coupled node
A1 = WE*e(m) - wi(1)*i(1); %note that we need to use the mth component here.  
A2 = theta*e(1); 
dy(1) = (-e(1) + 1./(1+exp(-a*(A1))))/tau1;
dy(m+1) = (-i(1) + 1./(1+exp(-a*A2)));
dy(2*m+1) = i(1)*(e(1)-p)/tau2;
dy(2:m) = D(2:m,:)*e; 
dy(m+2:2*m) = D(2:m,:)*i;
dy(2*m+2:3*m) = D(2:m,:)*wi;


%Linearized system
phiy = 1./(1+exp(-a*(WE*e(m)-wi(1)*i(1)))); 
phiPy = a*phiy*(1-phiy); 
phiY = 1./(1+exp(-a*(theta*e(1))));
phiPY = a*phiY*(1-phiY); 



J = zeros(3*m,3*m); 
J(1,1) = -1/tau1; 
J(1,m) = r*phiPy/tau1; 
J(1,m+1) = -wi(1)*phiPy/tau1 ;
J(1,2*m+1) = -i(1)*phiPy/tau1; 
J(2:m,1:m) = D(2:m,:); 
J(m+1,1) = theta*phiPY;
J(m+1,m+1) = -1; 
J(m+2:2*m,m+1:2*m) = D(2:m,:); 
J(2*m+1,1) = i(1)/tau2;
J(2*m+1,m+1) = (e(1)-p)/tau2;
J(2*m+2:3*m,2*m+1:3*m) = D(2:m,:);

dy(3*m+1:3*m+(3*m)^2) = J*Y;

end