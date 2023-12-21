function [dx, dy, dz] = CentralDerivative_d6_2ndOrder(u)


% Allocate memory for derivatives
dx = zeros(length(u(:,1,1)),length(u(1,:,1)),length(u(1,1,:)));
dy = dx;
dz = dx;

%% X direction
% Extremes left to 0 (Boundary points not needed)
% Sweep internal points [1 -6 15 -20
for i = 4:(length(u(1,:,1))-3)
    dx(:,i,:)    = (u(:,i+3,:) - 6*u(:,i+2,:) + 15*u(:,i+1,:) - 20*u(:,i,:) + 15*u(:,i-1,:) - 6*u(:,i-2,:)  + u(:,i-3,:));
end

% Forward / backward extremes
i = 2;
dx(:,i,:)    = (u(:,i,:) - 6*u(:,i+1,:) + 15*u(:,i+2,:) - 20*u(:,i+3,:) + 15*u(:,i+4,:) - 6*u(:,i+5,:)  + u(:,i+6,:));
i = 3;
dx(:,i,:)    = (u(:,i,:) - 6*u(:,i+1,:) + 15*u(:,i+2,:) - 20*u(:,i+3,:) + 15*u(:,i+4,:) - 6*u(:,i+5,:)  + u(:,i+6,:));

i = length(u(1,:,1))-2;
dx(:,i,:)    = (u(:,i,:) - 6*u(:,i-1,:) + 15*u(:,i-2,:) - 20*u(:,i-3,:) + 15*u(:,i-4,:) - 6*u(:,i-5,:)  + u(:,i-6,:));
i = length(u(1,:,1))-1;
dx(:,i,:)    = (u(:,i,:) - 6*u(:,i-1,:) + 15*u(:,i-2,:) - 20*u(:,i-3,:) + 15*u(:,i-4,:) - 6*u(:,i-5,:)  + u(:,i-6,:));

%% Y direction
for j = 4:(length(u(:,1,1))-3)
    dy(j,:,:)    = (u(j+3,:,:) - 6*u(j+2,:,:) + 15*u(j+1,:,:) - 20*u(j,:,:) + 15*u(j-1,:,:) - 6*u(j-2,:,:)  + u(j-3,:,:));
end

% Forward / backward extremes
j = 2;
dy(j,:,:)    = (u(j,:,:) - 6*u(j+1,:,:) + 15*u(j+2,:,:) - 20*u(j+3,:,:) + 15*u(j+4,:,:) - 6*u(j+5,:,:)  + u(j+6,:,:));
j = 3;
dy(j,:,:)    = (u(j,:,:) - 6*u(j+1,:,:) + 15*u(j+2,:,:) - 20*u(j+3,:,:) + 15*u(j+4,:,:) - 6*u(j+5,:,:)  + u(j+6,:,:));

j = length(u(:,1,1))-2;
dy(j,:,:)    = (u(j,:,:) - 6*u(j-1,:,:) + 15*u(j-2,:,:) - 20*u(j-3,:,:) + 15*u(j-4,:,:) - 6*u(j-5,:,:)  + u(j-6,:,:));
j = length(u(:,1,1))-1;
dy(j,:,:)    = (u(j,:,:) - 6*u(j-1,:,:) + 15*u(j-2,:,:) - 20*u(j-3,:,:) + 15*u(j-4,:,:) - 6*u(j-5,:,:)  + u(j-6,:,:));

%% Z direction
for k = 4:(length(u(1,1,:))-3)
    dz(:,:,k)    =  (u(:,:,k+3) - 6*u(:,:,k+2) + 15*u(:,:,k+1) - 20*u(:,:,k) + 15*u(:,:,k-1) - 6*u(:,:,k-2)  + u(:,:,k-3));
end

% Forward / backward extremes
k = 2;
dz(:,:,k)    =  (u(:,:,k) - 6*u(:,:,k+1) + 15*u(:,:,k+2) - 20*u(:,:,k+3) + 15*u(:,:,k+4) - 6*u(:,:,k+5)  + u(:,:,k+6));
k = 3;
dz(:,:,k)    =  (u(:,:,k) - 6*u(:,:,k+1) + 15*u(:,:,k+2) - 20*u(:,:,k+3) + 15*u(:,:,k+4) - 6*u(:,:,k+5)  + u(:,:,k+6));

k = length(u(1,1,:))-2;
dz(:,:,k)    =  (u(:,:,k) - 6*u(:,:,k-1) + 15*u(:,:,k-2) - 20*u(:,:,k-3) + 15*u(:,:,k-4) - 6*u(:,:,k-5)  + u(:,:,k-6));
k = length(u(1,1,:))-1;
dz(:,:,k)    =  (u(:,:,k) - 6*u(:,:,k-1) + 15*u(:,:,k-2) - 20*u(:,:,k-3) + 15*u(:,:,k-4) - 6*u(:,:,k-5)  + u(:,:,k-6));




end