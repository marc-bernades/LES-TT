function value_filt = FilterFields_1D_CDLF_Box(value,Epsilon_current,dx)

value_filt   = zeros(size(dx));
V_o          = zeros(size(dx));

% 1D in X direction (Compressible HPFC notation)
i = 2; k = 2;

try
    % Reference Baez Vidal et al. 2016 JCP
    for j = 2:length(value(1,:,1))-1

        % Adjust order near the borders
        if i <= Epsilon_current/2
            idx = 1;
        elseif i >= length(value(1,:,1)) - Epsilon_current/2
            idx = 1;
        else
            % Ratio doble than current epsilon
            idx = Epsilon_current/2;
        end

        V_o(i,j,k)         = dx(i,j,k);

        % i+1, j, k
        volume     = dx(i,j+idx,k);
        f_op(1)    = value(i,j+idx,k)*volume;
        f_op2(1)   = volume;

        % i-1, j, k
        volume     = dx(i,j-idx,k);
        f_op(2)    = value(i,j-idx,k)*volume;
        f_op2(2)   = volume;

    

        % Filt value
        value_filt(i,j,k) = 1/Epsilon_current*(value(i,j,k) + (Epsilon_current - 1)*sum(f_op)/sum(f_op2));

        % isnan check
        if isnan(value_filt(i,j,k))
            value_filt(i,j,k) = 0;
        end

    end

catch ME
    disp(ME.message)
end

end