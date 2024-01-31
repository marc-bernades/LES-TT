function [wn,ke_mag] = Calculate_TKE_spectra_3D(u,v,w,Lx,Ly,Lz,Nx,Ny,Nz)

% Grid in physical and Fourier space
N = (Nx*Ny*Nz)^(1/3);
m = -N/2:N/2-1;

% Reshape data in case it has changed
u = reshape(u,int16([N,N,N]));
v = reshape(v,int16([N,N,N]));
w = reshape(w,int16([N,N,N]));

% Wavenumber
kx = m*2*pi/Lx;
ky = m*2*pi/Ly;
kz = m*2*pi/Lz;
k  = sqrt(kx.^2 + ky.^2 + kz.^2);

% Fourier and shift wave numbers
U = fftn(u)/N^3; U = fftshift(U);
V = fftn(v)/N^3; V = fftshift(V);
W = fftn(w)/N^3; W = fftshift(W);


% Compute the kinetic energy
ke = 1/2*(U.*conj(U) + V.*conj(V) + W.*conj(W));

% Make the energy spectrum
m_x = length(ke(:,1,1)); m_y = length(ke(1,:,1)); m_z = length(ke(1,1,:));
m_radius = ceil(1/2*sqrt(m_x^2 + m_y^2 + m_z^2)) + 1;
epsilon  = 1E-50; % Avoid log(0)
ke_mag   = zeros(1,m_radius) + epsilon;

% Sides
center_x = m_x/2;
center_y = m_y/2;
center_z = m_z/2;


for ii = 1:m_x
    for jj = 1:m_y
        for kk = 1:m_z
            % wn
            wn = 1 + round(sqrt((ii-center_x)^2 + (jj-center_y)^2 + (kk-center_z)^2));
            ke_mag(wn) = ke_mag(wn) + ke(ii,jj,kk);
        end
    end
end

wn = (1:m_radius);


end