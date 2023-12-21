function value_filt = FilterFields(value,delta_filt,delta,Epsilon_current,x,y,z,dx,dy,dz,Filter_type)

value_filt   = zeros(size(dx));
V_o          = zeros(size(dx));
Epsilon_LDE  = zeros(size(dx)) + Epsilon_current;

switch Filter_type

    case 'Top_hat'

        % Order Taylor expansion to Stencil

        % 2nd order Taylor expansion
        if Epsilon_current == 2

            [d2u_y,d2u_x,d2u_z]  = CentralDerivative_d2_2ndOrder(value); % Order on h5 is X,Y,Z hence need to swapp
            Lap_U                = d2u_x./dx.^2 + d2u_y./dy.^2 + d2u_z./dz.^2;
            value_filt(2:end-1,2:end-1,2:end-1)  = value(2:end-1,2:end-1,2:end-1) + (delta_filt(2:end-1,2:end-1,2:end-1).^2)/24.*Lap_U(2:end-1,2:end-1,2:end-1);

            % 4th order Taylor expansion
        elseif Epsilon_current == 4

            [d2u_y,d2u_x,d2u_z]  = CentralDerivative_d2_2ndOrder(value); % Order on h5 is X,Y,Z hence need to swapp
            Lap_U                = d2u_x./dx.^2 + d2u_y./dy.^2 + d2u_z./dz.^2;
            [d4u_y,d4u_x,d4u_z]  = CentralDerivative_d4_2ndOrder(value);
            D4_U                 = d4u_x./dx.^4 + d4u_y./dy.^4 + d4u_z./dz.^4;
            value_filt(2:end-1,2:end-1,2:end-1)  = value(2:end-1,2:end-1,2:end-1) + (delta_filt(2:end-1,2:end-1,2:end-1).^2)/24.*Lap_U(2:end-1,2:end-1,2:end-1) + ...
                (delta_filt(2:end-1,2:end-1,2:end-1).^4)/1920.*D4_U(2:end-1,2:end-1,2:end-1);

            % 6th order Taylor expansion
        elseif Epsilon_current >= 6

            [d2u_y,d2u_x,d2u_z]  = CentralDerivative_d2_2ndOrder(value); % Order on h5 is X,Y,Z hence need to swapp
            Lap_U                = d2u_x./dx.^2 + d2u_y./dy.^2 + d2u_z./dz.^2;
            [d4u_y,d4u_x,d4u_z]  = CentralDerivative_d4_2ndOrder(value);
            D4_U                 = d4u_x./dx.^4 + d4u_y./dy.^4 + d4u_z./dz.^4;
            [d6u_y,d6u_x,d6u_z]  = CentralDerivative_d6_2ndOrder(value);
            D6_U                 = d6u_x./dx.^6 + d6u_y./dy.^6 + d6u_z./dz.^6;

            value_filt(2:end-1,2:end-1,2:end-1)  = value(2:end-1,2:end-1,2:end-1) + (delta_filt(2:end-1,2:end-1,2:end-1).^2)/24.*Lap_U(2:end-1,2:end-1,2:end-1) + ...
                (delta_filt(2:end-1,2:end-1,2:end-1).^4)/1920.*D4_U(2:end-1,2:end-1,2:end-1) + (delta_filt(2:end-1,2:end-1,2:end-1).^6)/322560.*D6_U(2:end-1,2:end-1,2:end-1);

        end



    case 'Top_hat_2_Forward'

        [d2u_y,d2u_x,d2u_z]  = ForwardDerivative_d2_2ndOrder(value); % Order on h5 is X,Y,Z hence need to swapp
        Lap_U                = d2u_x./dx.^2 + d2u_y./dy.^2 + d2u_z./dz.^2;
        value_filt(2:end-1,2:end-1,2:end-1)  = value(2:end-1,2:end-1,2:end-1) + (delta_filt(2:end-1,2:end-1,2:end-1).^2)/24.*Lap_U(2:end-1,2:end-1,2:end-1);

    case 'Digital_image'

        % Average in the 9x3 neightbor
        for i = 2:length(value(:,1,1))-1
            for j = 2:length(value(1,:,1))-1
                for k = 2:length(value(:,:,1))-1
                    neightbor(1)  = value(i+1,j,k);
                    neightbor(2)  = value(i-1,j,k);
                    neightbor(3)  = value(i,j+1,k);
                    neightbor(4)  = value(i,j-1,k);
                    neightbor(5)  = value(i,j,k+1);
                    neightbor(6)  = value(i,j,k-1);
                    neightbor(7)  = value(i+1,j-1,k);
                    neightbor(8)  = value(i-1,j-1,k);
                    neightbor(9)  = value(i,j-1,k+1);
                    neightbor(10) = value(i,j+1,k-1);
                    neightbor(11) = value(i+1,j+1,k);
                    neightbor(12) = value(i-1,j+1,k);
                    neightbor(13) = value(i,j+1,k+1);
                    neightbor(14) = value(i,j+1,k-1);
                    neightbor(15) = value(i+1,j,k+1);
                    neightbor(16) = value(i-1,j,k+1);
                    neightbor(17) = value(i,j+1,k+1);
                    neightbor(18) = value(i+1,j,k-1);
                    neightbor(19) = value(i-1,j,k-1);
                    neightbor(20) = value(i,j-1,k-1);
                    neightbor(21) = value(i+1,j+1,k+1);
                    neightbor(22) = value(i-1,j-1,k-1);
                    neightbor(23) = value(i+1,j-1,k+1);
                    neightbor(24) = value(i+1,j-1,k-1);
                    neightbor(25) = value(i-1,j+1,k+1);
                    neightbor(26) = value(i-1,j-1,k+1);
                    neightbor(27) = value(i,j,k);
                    value_filt(i,j,k) = sum(neightbor)/length(neightbor);


                end
            end
        end


    case 'CDLF_Lap'

        % Reference Baez Vidal et al. 2016 JCP
        for i = 2:length(value(:,1,1))-1
            for j = 2:length(value(1,:,1))-1
                for k = 2:length(value(:,:,1))-1
                    V_o(i,j,k)         = dx(i,j,k)*dy(i,j,k)*dz(i,j,k);

                    % i+1, j, k
                    A_op       = dy(i,j,k)*dz(i,j,k);
                    n_op       = [1,0,0];
                    r_op       = [x(i+1,j,k) - x(i,j,k),y(i+1,j,k) - y(i,j,k),z(i+1,j,k) - z(i,j,k)];
                    f_op(1)    = (value(i+1,j,k) - value(i,j,k))*A_op./(dot(n_op,r_op));
                    f_op2(1)   = A_op./(dot(n_op,r_op));

                    % i-1, j, k
                    A_op       = dy(i,j,k)*dz(i,j,k);
                    n_op       = [-1,0,0];
                    r_op       = [x(i-1,j,k) - x(i,j,k),y(i-1,j,k) - y(i,j,k),z(i-1,j,k) - z(i,j,k)];
                    f_op(2)    = (value(i-1,j,k) - value(i,j,k))*A_op./(dot(n_op,r_op));
                    f_op2(2)   = A_op./(dot(n_op,r_op));

                    % i, j+1, k
                    A_op       = dx(i,j,k)*dz(i,j,k);
                    n_op       = [0,1,0];
                    r_op       = [x(i,j+1,k) - x(i,j,k),y(i,j+1,k) - y(i,j,k),z(i,j+1,k) - z(i,j,k)];
                    f_op(3)    = (value(i,j+1,k) - value(i,j,k))*A_op./(dot(n_op,r_op));
                    f_op2(3)   = A_op./(dot(n_op,r_op));

                    % i, j-1, k
                    A_op       = dx(i,j,k)*dz(i,j,k);
                    n_op       = [0,-1,0];
                    r_op       = [x(i,j-1,k) - x(i,j,k),y(i,j-1,k) - y(i,j,k),z(i,j-1,k) - z(i,j,k)];
                    f_op(4)    = (value(i,j-1,k) - value(i,j,k))*A_op./(dot(n_op,r_op));
                    f_op2(4)   = A_op./(dot(n_op,r_op));

                    % i, j, k+1
                    A_op       = dx(i,j,k)*dy(i,j,k);
                    n_op       = [0,0,1];
                    r_op       = [x(i,j,k+1) - x(i,j,k),y(i,j,k+1) - y(i,j,k),z(i,j,k+1) - z(i,j,k)];
                    f_op(5)    = (value(i,j,k+1) - value(i,j,k))*A_op./(dot(n_op,r_op));
                    f_op2(5)   = A_op./(dot(n_op,r_op));

                    % i, j, k-1
                    A_op       = dx(i,j,k)*dy(i,j,k);
                    n_op       = [0,0,-1];
                    r_op       = [x(i,j,k-1) - x(i,j,k),y(i,j,k-1) - y(i,j,k),z(i,j,k-1) - z(i,j,k)];
                    f_op(6)    = (value(i,j,k-1) - value(i,j,k))*A_op./(dot(n_op,r_op));
                    f_op2(6)   = A_op./(dot(n_op,r_op));

                    % Criteria
                    Sigma    = 1;
                    LED_Cond = 24*V_o(i,j,k)^(1/3)/(Sigma*sum(f_op2));
                    if Epsilon_current^2 > LED_Cond
                        Epsilon_LDE(i,j,k) = sqrt(LED_Cond);
                    end


                    % Filt value
                    value_filt(i,j,k) = value(i,j,k) + (Epsilon_LDE(i,j,k)*delta(i,j,k)).^2/(24*V_o(i,j,k))*sum(f_op);

                end
            end
        end

    case 'CDLF_Box'

        try
        % Reference Baez Vidal et al. 2016 JCP
        for i = 2:length(value(:,1,1))-1
            for j = 2:length(value(1,:,1))-1
                for k = 2:length(value(1,1,:))-1 %length(value(:,:,1))-1
 
                    % Adjust order near the borders
                    if i <= Epsilon_current/2 || j <= Epsilon_current/2 || k <= Epsilon_current/2
                        idx = 1;
                    elseif i >= length(value(:,1,1)) - Epsilon_current/2 || j >= length(value(1,:,1)) - Epsilon_current/2 || k >= length(value(1,1,:)) - Epsilon_current/2
                        idx = 1;
                    else
                        % Ratio doble than current epsilon
                        idx = Epsilon_current/2;
                    end

                    V_o(i,j,k)         = dx(i,j,k)*dy(i,j,k)*dz(i,j,k);

                    % i+1, j, k
                    volume     = dx(i+idx,j,k)*dy(i+idx,j,k)*dz(i+idx,j,k);
                    f_op(1)    = value(i+idx,j,k)*volume;
                    f_op2(1)   = volume;

                    % i-1, j, k
                    volume     = dx(i-idx,j,k)*dy(i-idx,j,k)*dz(i-idx,j,k);
                    f_op(2)    = value(i-idx,j,k)*volume;
                    f_op2(2)   = volume;

                    % i, j+1, k
                    volume     = dx(i,j+idx,k)*dy(i,j+idx,k)*dz(i,j+idx,k);
                    f_op(3)    = value(i,j+idx,k)*volume;
                    f_op2(3)   = volume;

                    % i, j-1, k
                    volume     = dx(i,j-idx,k)*dy(i,j-idx,k)*dz(i,j-idx,k);
                    f_op(4)    = value(i,j-idx,k)*volume;
                    f_op2(4)   = volume;

                    % i, j, k+1
                    volume     = dx(i,j,k+idx)*dy(i,j,k+idx)*dz(i,j,k+idx);
                    f_op(5)    = value(i,j,k+idx)*volume;
                    f_op2(5)   = volume;

                    % i, j, k-1
                    volume     = dx(i,j,k-idx)*dy(i,j,k-idx)*dz(i,j,k-idx);
                    f_op(6)    = value(i,j,k-idx)*volume;
                    f_op2(6)   = volume;

                    % Filt value
                    value_filt(i,j,k) = 1/Epsilon_current*(value(i,j,k) + (Epsilon_current - 1)*sum(f_op)/sum(f_op2));

                    % isnan check
                    if isnan(value_filt(i,j,k))
                        value_filt(i,j,k) = 0;
                    end

                end
            end
        end

        catch ME
            disp(ME.message)
        end

end