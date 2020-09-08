function [Sum_I,Sum_Q,CodeNcoPhase,CarrierNcoPhase,CodeIdx] = CorherentSum(Input_I,Input_Q,CAcode1023,CarrierNcoStep,CodeNcoStep,CarrierNcoPhase...
    ,CodeNcoPhase,Ms,sum_I,sum_Q,codeidx)
%CORHERENTSUM 此处显示有关此函数的摘要
%   此处显示详细说明
global sin_table;
global cos_table;

    Sum_I = sum_I;
    Sum_Q = sum_Q;
    CodeIdx = codeidx;
    CodeNcoPhaseOld = CodeNcoPhase;

    for n=1:1:2048*Ms
        CodeNcoPhase = mod(CodeNcoPhase + CodeNcoStep,2^32-1);
        CarrierNcoPhase = mod(CarrierNcoPhase + CarrierNcoStep,2^32-1);
        PhaseIdx = bitand(bitshift(uint32(CarrierNcoPhase),-28),15)+1;
        I_temp = Input_I(n) * cos_table(PhaseIdx) - Input_Q(n) * sin_table(PhaseIdx);
        Q_temp = Input_Q(n) * cos_table(PhaseIdx) + Input_I(n) * sin_table(PhaseIdx);
        Sum_I = Sum_I + I_temp * CAcode1023(mod(CodeIdx,1023)+1);
        Sum_Q = Sum_Q + Q_temp * CAcode1023(mod(CodeIdx,1023)+1);

        if(CodeNcoPhase<CodeNcoPhaseOld)
            CodeIdx = CodeIdx + 1;
        end
        CodeNcoPhaseOld = CodeNcoPhase;

    end
end

