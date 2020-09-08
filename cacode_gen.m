function [CA_CODE] = cacode_gen(prn);
    CA_CODE = zeros(1,1023);
    G1_reg = ones(1,10);

    % G2_reg = [0 0 0 1 0 1 0 0 0 1];% QZSS 193 MSB -> LSB
    % G2_reg = [1 1 1 0 0 0 0 1 1 1];% QZSS 194
    table = [0 0 0 0 0 1 0 0 1 1
    0 0 0 0 1 0 0 1 1 1
    0 0 0 1 0 0 1 1 1 1
    0 0 1 0 0 1 1 1 1 1
    1 1 0 1 1 0 1 0 0 1
    1 0 1 1 0 1 0 0 1 1
    1 0 0 1 1 0 1 0 0 1
    0 0 1 1 0 1 0 0 1 1
    0 1 1 0 1 0 0 1 1 1
    0 0 1 0 0 0 1 0 1 1
    0 1 0 0 0 1 0 1 1 1
    0 0 0 1 0 1 1 1 1 1
    0 0 1 0 1 1 1 1 1 1
    0 1 0 1 1 1 1 1 1 1
    1 0 1 1 1 1 1 1 1 1
    0 1 1 1 1 1 1 1 1 1 %GPS16
    0 1 1 1 0 1 1 0 0 1
    1 1 1 0 1 1 0 0 1 1
    1 1 0 1 1 0 0 1 1 1
    1 0 1 1 0 0 1 1 1 1
    0 1 1 0 0 1 1 1 1 1
    1 1 0 0 1 1 1 1 1 1];%GPS22
    G2_reg = ~table(prn,:);
    for i=1:1:1023,
        G1_tmp_reg = mod(G1_reg(3) + G1_reg(10),2);
        G2_tmp_reg = mod(G2_reg(2) + G2_reg(3) + G2_reg(6) + G2_reg(8) + G2_reg(9)...
           + G2_reg(10),2);
        CA_CODE(i) = mod(G1_reg(10) + G2_reg(10),2);
        G1_reg = [G1_tmp_reg G1_reg(1:9)];
        G2_reg = [G2_tmp_reg G2_reg(1:9)];
    end
    CA_CODE = (2*CA_CODE - 1);
end
