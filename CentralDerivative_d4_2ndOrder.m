function [dx, dy, dz] = CentralDerivative_d4_2ndOrder(u)


% Allocate memory for derivatives
dx = zeros(length(u(:,1,1)),length(u(1,:,1)),length(u(1,1,:)));
dy = dx;
dz = dx;

%% X direction
% Extremes left to 0 (Boundary points not needed)
% Sweep internal points
for i = 3:(length(u(1,:,1))-2)
    dx(:,i,:)    = (u(:,i+2,:) - 4*u(:,i+1,:) + 6*u(:,i,:) - 4*u(:,i-1,:) + u(:,i-2,:));
end

% Forward / backward extremes
i = 2;
dx(:,i,:)    = (u(:,i,:) - 4*u(:,i+1,:) + 6*u(:,i+2,:) - 4*u(:,i+3,:) + u(:,i+4,:));
i = length(u(1,:,1))-1;
dx(:,i,:)    = (u(:,i,:) - 4*u(:,i-1,:) + 6*u(:,i-2,:) - 4*u(:,i-3,:) + u(:,i-4,:));

%% Y direction
for j = 3:(length(u(:,1,1))-2)
    dy(j,:,:)    = (u(j+2,:,:) - 4*u(j+1,:,:) + 6*u(j,:,:) - 4*u(j-1,:,:) + u(j-2,:,:));
end

% Forward / backward extremes
j = 2;
dy(j,:,:)    = (u(j,:,:) - 4*u(j+1,:,:) + 6*u(j+2,:,:) - 4*u(j+3,:,:) + u(j+4,:,:));
j = length(u(:,1,1))-1;
dy(j,:,:)    = (u(j,:,:) - 4*u(j-1,:,:) + 6*u(j-2,:,:) - 4*u(j-3,:,:) + u(j-4,:,:));


%% Z direction
for k = 3:(length(u(1,1,:))-2)
    dz(:,:,k)    = (u(:,:,k+2) - 4*u(:,:,k+1) + 6*u(:,:,k)  - 4*u(:,:,k-1) + u(:,:,k-2));
end

% Forward / backward extremes
k = 2;
dz(:,:,k)    = (u(:,:,k) - 4*u(:,:,k+1) + 6*u(:,:,k+2)  - 4*u(:,:,k+3) + u(:,:,k+4));
k = length(u(1,1,:))-1;
dz(:,:,k)    = (u(:,:,k) - 4*u(:,:,k-1) + 6*u(:,:,k-2)  - 4*u(:,:,k-3) + u(:,:,k-4));



end