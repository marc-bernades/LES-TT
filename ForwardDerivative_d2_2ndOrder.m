function [dx, dy, dz] = ForwardDerivative_d2_2ndOrder(u)


% Allocate memory for derivatives
dx = zeros(length(u(:,1,1)),length(u(1,:,1)),length(u(1,1,:)));
dy = dx;
dz = dx;

%% X direction
% Extremes left to 0 (Boundary points not needed)
% Sweep internal points
for i = 2:(length(u(1,:,1))-2)
    dx(:,i,:)    = (u(:,i+2,:) - 2*u(:,i+1,:) + u(:,i,:));
end

% Central last point - 1
i = length(u(1,:,1))-1;
dx(:,i,:)    = (u(:,i+1,:) - 2*u(:,i,:) + u(:,i-1,:));


%% Y direction
for j = 2:(length(u(:,1,1))-2)
    dy(j,:,:)    = (u(j+2,:,:) - 2*u(j+1,:,:) + u(j,:,:));
end

% Central last point - 1
j = length(u(:,1,1))-1;
dy(j,:,:)    = (u(j+1,:,:) - 2*u(j,:,:) + u(j-1,:,:));


%% Z direction
for k = 2:(length(u(1,1,:))-2)
    dz(:,:,k)    = (u(:,:,k+2) - 2*u(:,:,k+1) + u(:,:,k));
end

% Central last point - 1
k = length(u(1,1,:))-1;
dz(:,:,k)    = (u(:,:,k+1) - 2*u(:,:,k) + u(:,:,k-1));



end

