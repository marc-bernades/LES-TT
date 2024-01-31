function [k_mag,e_k_mag] = Calculate_TKE_spectra(u,v,w,L,N)

% N=num_points_z;
% delta = 100*1E-6;
% Lx=4*delta; %length in the z direction
% Lz=4/3*delta; %length in the z direction
% 
% kz    = 2*pi*(-N/2:1:N/2-1)/Lz; %frequencies, or wavenumber
% Uzre  = w(:,2,:);%change the column of Ux into a matrix,each column corresponding to one z
% Uzmre = avg_w(:,2,:);
% Uzmreave=mean(Uzmre(:));%field average
% Uzmfield=zeros(N,N);
%  for i=1:N
%      for j=1:N
%          Uzmfield(i,j)=Uzmreave;
%      end
%  end
%  w_k=Uzre-Uzmfield;%get the velocity fluctuation=instannous velocity - timeaveraged velocity
%  wk=fft(w_k);%each column corresponds to a z value
%  wkk=fftshift(wk);
%  E_k=zeros(N,N);
% for i=1:N
%     for j=1:N
%     E_k(i,j)=(wkk(i,j)*conj(wkk(i,j)));%equal to abs(wkmean(i))^2
%     end
% end
% Emean=mean(E_k,2);
% loglog(kz,Emean);


% Make grid s in physical and Fourier space
m = -N/2:N/2-1;
m = m*2*pi/L;

% Fourier and shift wave numbers
% perform the Fourier transform
U_hat = fftn(u);
V_hat = fftn(v);
W_hat = fftn(w);

% shift the Fourier transform
U_hat = fftshift(U_hat);
V_hat = fftshift(V_hat);
W_hat = fftshift(W_hat);

%% Generate a spectrum of the Fourier transforms with a given resolution

% make the physical Fourier transforms
U_hat_phys = (L/N)^3*U_hat;
V_hat_phys = (L/N)^3*V_hat;
W_hat_phys = (L/N)^3*W_hat;

% Compute the energy of the field
e_phys = 1/2*(U_hat_phys.*conj(U_hat_phys) + V_hat_phys.*conj(V_hat_phys) + W_hat_phys.*conj(W_hat_phys));

% Make the energy spectrum
m_max   = ceil(sqrt(3)*N/2);
e_k_mag = zeros(1,m_max+1);


for im = 1:N
    for jm = 1:N
        for km = 1:N
            % make a grid of |m| values around each m = (im,jm,km)
            mgrid = reshape(round(sqrt( (m(im)).^2 + (m(jm)).^2 + (m(km)).^2 )),1,1);
            % add the energy for each |m|
            e_k_per_cell = e_phys(im,jm,km);
            e_k_mag(mgrid+1) = e_k_mag(mgrid+1) + e_k_per_cell;
        end
    end
end

% compute the energy per physical wave number k (instead of m)
k_mag = 2*pi/L*(0:m_max);
e_k_mag = e_k_mag*L/(2*pi);


% dk = 2*pi/L;
% Ef = sum(e_k_mag)*dk;
% epsf = sum(2*nu.*(k_mag/dk).^2.*e_k_mag)*dk;
% Rl = sqrt(20/3)*Ef*1/sqrt(nu*epsf);
% kol = (nu.^3/epsf)^(1/4);
% ensc = (epsf*nu^5)^(1/4);

k_mag   = k_mag/2/pi;
e_k_mag = e_k_mag*2*pi;




end