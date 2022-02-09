#ifndef WINSLOWCSEARRAY_H_IS_INCLUDED
#define WINSLOWCSEARRAY_H_IS_INCLUDED

#include <memory>
#include <stdexcept>
#include <cmath>

#include "goss/ParameterizedODE.h"

namespace goss {

  // Implementation of gotran generated ODE
  class WinslowCSEArray : public ParameterizedODE
  {
  public:

    // Constructor
    WinslowCSEArray() : ParameterizedODE(31, 68, 1, 2, 4),
      A_cap(0.0001534), C_sc(1.0), V_JSR(1.6e-07), V_NSR(2.1e-06),
        V_myo(2.584e-05), V_ss(1.2e-09), ist(0), Ca_o(2.0), K_o(4.0),
        Na_i(10.0), Na_o(138.0), G_KpMax(0.002216), G_KrMax(0.0034),
        G_KsMax(0.00271), G_NaMax(12.8), G_bCaMax(0.0003842),
        G_bNaMax(0.0031), G_tiMax(2.8), G_toMax(0.23815), I_NaKMax(0.693),
        I_pCaMax(0.05), K_mCa(1.38), K_mK1(13.0), K_mKo(1.5), K_mNa(87.5),
        K_mNai(10.0), K_mpCa(5e-05), eta(0.35), k_NaCa(0.3), k_sat(0.2),
        K_SR(1.0), K_fb(0.000168), K_rb(3.29), N_fb(1.2), N_rb(1.0),
        kaminus(0.576), kaplus(0.01215), kbminus(1.93), kbplus(0.00405),
        kcminus(0.0008), kcplus(0.1), mcoop(3.0), ncoop(4.0), tau_tr(0.5747),
        tau_xfer(26.7), v_1(1.8), v_maxf(8.13e-05), v_maxr(0.000318),
        ICahalf(-0.265), PCa(0.0003125), PK(5.79e-07), aL(2.0), bL(2.0),
        fL(0.3), gL(2.0), omega(0.01), CMDNtot(0.05), CSQNtot(15.0),
        EGTAtot(10.0), HTRPNtot(0.14), KmCMDN(0.00238), KmCSQN(0.8),
        KmEGTA(0.00015), LTRPNtot(0.07), khtrpn_minus(6.6e-05),
        khtrpn_plus(20.0), kltrpn_minus(0.04), kltrpn_plus(40.0)

    {

      // State names
      _state_names[0] = "V";
      _state_names[1] = "h";
      _state_names[2] = "j";
      _state_names[3] = "m";
      _state_names[4] = "xKr";
      _state_names[5] = "xKs";
      _state_names[6] = "xto1";
      _state_names[7] = "yto1";
      _state_names[8] = "K_i";
      _state_names[9] = "Ca_JSR";
      _state_names[10] = "Ca_NSR";
      _state_names[11] = "Ca_i";
      _state_names[12] = "Ca_ss";
      _state_names[13] = "C1_RyR";
      _state_names[14] = "C2_RyR";
      _state_names[15] = "O1_RyR";
      _state_names[16] = "O2_RyR";
      _state_names[17] = "C0";
      _state_names[18] = "C1";
      _state_names[19] = "C2";
      _state_names[20] = "C3";
      _state_names[21] = "C4";
      _state_names[22] = "CCa0";
      _state_names[23] = "CCa1";
      _state_names[24] = "CCa2";
      _state_names[25] = "CCa3";
      _state_names[26] = "CCa4";
      _state_names[27] = "Open";
      _state_names[28] = "yCa";
      _state_names[29] = "HTRPNCa";
      _state_names[30] = "LTRPNCa";

      // Parameter names
      _parameter_names[0] = "A_cap";
      _parameter_names[1] = "C_sc";
      _parameter_names[2] = "V_JSR";
      _parameter_names[3] = "V_NSR";
      _parameter_names[4] = "V_myo";
      _parameter_names[5] = "V_ss";
      _parameter_names[6] = "ist";
      _parameter_names[7] = "Ca_o";
      _parameter_names[8] = "K_o";
      _parameter_names[9] = "Na_i";
      _parameter_names[10] = "Na_o";
      _parameter_names[11] = "G_KpMax";
      _parameter_names[12] = "G_KrMax";
      _parameter_names[13] = "G_KsMax";
      _parameter_names[14] = "G_NaMax";
      _parameter_names[15] = "G_bCaMax";
      _parameter_names[16] = "G_bNaMax";
      _parameter_names[17] = "G_tiMax";
      _parameter_names[18] = "G_toMax";
      _parameter_names[19] = "I_NaKMax";
      _parameter_names[20] = "I_pCaMax";
      _parameter_names[21] = "K_mCa";
      _parameter_names[22] = "K_mK1";
      _parameter_names[23] = "K_mKo";
      _parameter_names[24] = "K_mNa";
      _parameter_names[25] = "K_mNai";
      _parameter_names[26] = "K_mpCa";
      _parameter_names[27] = "eta";
      _parameter_names[28] = "k_NaCa";
      _parameter_names[29] = "k_sat";
      _parameter_names[30] = "K_SR";
      _parameter_names[31] = "K_fb";
      _parameter_names[32] = "K_rb";
      _parameter_names[33] = "N_fb";
      _parameter_names[34] = "N_rb";
      _parameter_names[35] = "kaminus";
      _parameter_names[36] = "kaplus";
      _parameter_names[37] = "kbminus";
      _parameter_names[38] = "kbplus";
      _parameter_names[39] = "kcminus";
      _parameter_names[40] = "kcplus";
      _parameter_names[41] = "mcoop";
      _parameter_names[42] = "ncoop";
      _parameter_names[43] = "tau_tr";
      _parameter_names[44] = "tau_xfer";
      _parameter_names[45] = "v_1";
      _parameter_names[46] = "v_maxf";
      _parameter_names[47] = "v_maxr";
      _parameter_names[48] = "ICahalf";
      _parameter_names[49] = "PCa";
      _parameter_names[50] = "PK";
      _parameter_names[51] = "aL";
      _parameter_names[52] = "bL";
      _parameter_names[53] = "fL";
      _parameter_names[54] = "gL";
      _parameter_names[55] = "omega";
      _parameter_names[56] = "CMDNtot";
      _parameter_names[57] = "CSQNtot";
      _parameter_names[58] = "EGTAtot";
      _parameter_names[59] = "HTRPNtot";
      _parameter_names[60] = "KmCMDN";
      _parameter_names[61] = "KmCSQN";
      _parameter_names[62] = "KmEGTA";
      _parameter_names[63] = "LTRPNtot";
      _parameter_names[64] = "khtrpn_minus";
      _parameter_names[65] = "khtrpn_plus";
      _parameter_names[66] = "kltrpn_minus";
      _parameter_names[67] = "kltrpn_plus";

      // Field state names
      _field_state_names[0] = "V";

      // Field state indices
      _field_state_indices[0] = 0;

      // Field parameter names
      _field_parameter_names[0] = "G_KrMax";
      _field_parameter_names[1] = "G_KsMax";

      // Monitored names
      _monitored_names[0] = "I_Ca";
      _monitored_names[1] = "I_CaK";
      _monitored_names[2] = "I_NaCa";
      _monitored_names[3] = "I_pCa";

      // Parameter to value map
      _param_to_value["A_cap"] = &A_cap;
      _param_to_value["C_sc"] = &C_sc;
      _param_to_value["V_JSR"] = &V_JSR;
      _param_to_value["V_NSR"] = &V_NSR;
      _param_to_value["V_myo"] = &V_myo;
      _param_to_value["V_ss"] = &V_ss;
      _param_to_value["ist"] = &ist;
      _param_to_value["Ca_o"] = &Ca_o;
      _param_to_value["K_o"] = &K_o;
      _param_to_value["Na_i"] = &Na_i;
      _param_to_value["Na_o"] = &Na_o;
      _param_to_value["G_KpMax"] = &G_KpMax;
      _param_to_value["G_KrMax"] = &G_KrMax;
      _param_to_value["G_KsMax"] = &G_KsMax;
      _param_to_value["G_NaMax"] = &G_NaMax;
      _param_to_value["G_bCaMax"] = &G_bCaMax;
      _param_to_value["G_bNaMax"] = &G_bNaMax;
      _param_to_value["G_tiMax"] = &G_tiMax;
      _param_to_value["G_toMax"] = &G_toMax;
      _param_to_value["I_NaKMax"] = &I_NaKMax;
      _param_to_value["I_pCaMax"] = &I_pCaMax;
      _param_to_value["K_mCa"] = &K_mCa;
      _param_to_value["K_mK1"] = &K_mK1;
      _param_to_value["K_mKo"] = &K_mKo;
      _param_to_value["K_mNa"] = &K_mNa;
      _param_to_value["K_mNai"] = &K_mNai;
      _param_to_value["K_mpCa"] = &K_mpCa;
      _param_to_value["eta"] = &eta;
      _param_to_value["k_NaCa"] = &k_NaCa;
      _param_to_value["k_sat"] = &k_sat;
      _param_to_value["K_SR"] = &K_SR;
      _param_to_value["K_fb"] = &K_fb;
      _param_to_value["K_rb"] = &K_rb;
      _param_to_value["N_fb"] = &N_fb;
      _param_to_value["N_rb"] = &N_rb;
      _param_to_value["kaminus"] = &kaminus;
      _param_to_value["kaplus"] = &kaplus;
      _param_to_value["kbminus"] = &kbminus;
      _param_to_value["kbplus"] = &kbplus;
      _param_to_value["kcminus"] = &kcminus;
      _param_to_value["kcplus"] = &kcplus;
      _param_to_value["mcoop"] = &mcoop;
      _param_to_value["ncoop"] = &ncoop;
      _param_to_value["tau_tr"] = &tau_tr;
      _param_to_value["tau_xfer"] = &tau_xfer;
      _param_to_value["v_1"] = &v_1;
      _param_to_value["v_maxf"] = &v_maxf;
      _param_to_value["v_maxr"] = &v_maxr;
      _param_to_value["ICahalf"] = &ICahalf;
      _param_to_value["PCa"] = &PCa;
      _param_to_value["PK"] = &PK;
      _param_to_value["aL"] = &aL;
      _param_to_value["bL"] = &bL;
      _param_to_value["fL"] = &fL;
      _param_to_value["gL"] = &gL;
      _param_to_value["omega"] = &omega;
      _param_to_value["CMDNtot"] = &CMDNtot;
      _param_to_value["CSQNtot"] = &CSQNtot;
      _param_to_value["EGTAtot"] = &EGTAtot;
      _param_to_value["HTRPNtot"] = &HTRPNtot;
      _param_to_value["KmCMDN"] = &KmCMDN;
      _param_to_value["KmCSQN"] = &KmCSQN;
      _param_to_value["KmEGTA"] = &KmEGTA;
      _param_to_value["LTRPNtot"] = &LTRPNtot;
      _param_to_value["khtrpn_minus"] = &khtrpn_minus;
      _param_to_value["khtrpn_plus"] = &khtrpn_plus;
      _param_to_value["kltrpn_minus"] = &kltrpn_minus;
      _param_to_value["kltrpn_plus"] = &kltrpn_plus;

    }

    // Evaluate rhs of the ODE
    void eval(const double* states, double t, double* values)
    {

      // Assign states

      // Common Sub Expressions
      const double cse_0 = K_o/states[8];
      const double cse_1 = (Na_o*Na_o*Na_o);
      const double cse_2 = 1.0/Na_i;
      const double cse_3 = 1.0/V_NSR;
      const double cse_4 = 1.0/tau_tr;
      const double cse_5 = states[10] - states[9];
      const double cse_6 = 1.0/V_myo;
      const double cse_7 = states[15] + states[16];
      const double cse_8 = -states[12] + states[9];
      const double cse_9 = 1.0/V_ss;
      const double cse_10 = 1.0/tau_xfer;
      const double cse_11 = -states[11] + states[12];
      const double cse_12 = 1000.0*states[12];
      const double cse_13 = 1.0/bL;
      const double cse_14 = (aL*aL);
      const double cse_15 = std::pow(cse_14, 3/2);
      const double cse_16 = std::pow(cse_15, 4/3);
      const double cse_17 = 0.2*states[0];
      const double cse_18 = 0.10375*states[12];
      const double cse_19 = A_cap*C_sc;
      const double cse_20 = kbminus*states[16];
      const double cse_21 = CMDNtot*KmCMDN;
      const double cse_22 = 0.0374417034617086*states[0];
      const double cse_23 = -states[7];
      const double cse_24 = -states[15];
      const double cse_25 = -states[29];
      const double cse_26 = -states[11];
      const double cse_27 = gL*states[27];
      const double cse_28 = cse_13*states[23];
      const double cse_29 = -states[1];
      const double cse_30 = -states[30];
      const double cse_31 = -states[2];
      const double cse_32 = kcminus*states[14];
      const double cse_33 = std::exp(0.0748834069234172*states[0]);
      const double cse_34 = std::exp(cse_22);
      const double cse_35 = -1.0*states[0];
      const double cse_36 = -0.1*states[0];
      const double cse_37 = 1.0/(K_mKo + K_o);
      const double cse_38 = std::pow(states[10]/K_rb, N_rb);
      const double cse_39 = std::pow(states[11]/K_fb, N_fb);
      const double cse_40 = 1.0/(K_mpCa + states[11]);
      const double cse_41 = 1.0/(Ca_o + K_mCa);
      const double cse_42 = std::pow(cse_12, ncoop);
      const double cse_43 = std::pow(cse_12, mcoop);
      const double cse_44 = (cse_13*cse_13);
      const double cse_45 = std::pow(cse_44, 3/2);
      const double cse_46 = std::pow(cse_45, 4/3);
      const double cse_47 = omega*cse_45;
      const double cse_48 = cse_10*cse_11;
      const double cse_49 = -1.0*K_o;
      const double cse_50 = cse_4*cse_5;
      const double cse_51 = omega*cse_44;
      const double cse_52 = omega*cse_46;
      const double cse_53 = I_pCaMax*cse_40;
      const double cse_54 = std::exp(0.1691*states[0] - 5.495);
      const double cse_55 = std::exp(cse_17 + 6.7);
      const double cse_56 = std::exp(0.1*states[0] + 0.2);
      const double cse_57 = kbplus*cse_43;
      const double cse_58 = states[0] - 26.7081865284974*std::log(cse_0);
      const double cse_59 = std::exp(cse_22*(eta - 1.0));
      const double cse_60 = 1.0/(cse_33 - 1.0);
      const double cse_61 = states[0] - 26.7081865284974*std::log(Na_o*cse_2);
      const double cse_62 = std::exp(-0.2*states[0] - 6.7);
      const double cse_63 = 1.0/(std::pow(K_mNai*cse_2, 1.5) + 1.0);
      const double cse_64 = states[0] -
        13.3540932642487*std::log(Ca_o/states[11]);
      const double cse_65 = 1.0/((K_mNa*K_mNa*K_mNa) + cse_1);
      const double cse_66 = std::exp(-0.0769230769230769*states[0] -
        0.153846153846154);
      const double cse_67 = 1.6*cse_56;
      const double cse_68 = aL*cse_56;
      const double cse_69 = 0.4*cse_56;
      const double cse_70 = 0.8*cse_56;
      const double cse_71 = 1.2*cse_56;
      const double cse_72 = kaplus*cse_42*states[13];
      const double cse_73 = v_1*cse_7*cse_8;
      const double cse_74 = -0.341*Ca_o + 0.001*cse_33;
      const double cse_75 = 1.0/(1.0 +
        0.1245*std::exp(-0.00374417034617086*states[0]));
      const double cse_76 = khtrpn_minus*cse_25 +
        khtrpn_plus*states[11]*(cse_25 + 1.0);
      const double cse_77 = kltrpn_minus*cse_30 +
        kltrpn_plus*states[11]*(cse_30 + 1.0);
      const double cse_78 = 0.1*cse_66;
      const double cse_79 = G_bCaMax*cse_64;
      const double cse_80 = 0.05*cse_66;
      const double cse_81 = 0.2*cse_66;
      const double cse_82 = 0.15*cse_66;
      const double cse_83 = 1.0/(std::exp(cse_35 - 40.0) + 1.0);
      const double cse_84 = 1.0/(k_sat*cse_59 + 1.0);
      const double cse_85 = 1.0/(cse_38 + cse_39 + 1.0);
      const double cse_86 = v_maxf*cse_39 - v_maxr*cse_38;
      const double cse_87 = G_toMax*cse_23*cse_58*states[6];
      const double cse_88 = 1.0/(cse_54 + std::exp(-0.0128*states[0] - 7.677));
      const double cse_89 = I_NaKMax*cse_37*cse_63*cse_75;
      const double cse_90 =
        -1.0*G_KpMax*cse_58/(std::exp(-0.167224080267559*states[0] +
        1.25217391304348) + 1.0);
      const double cse_91 = -G_KsMax*(states[5]*states[5])*(states[0] -
        26.7081865284974*std::log((K_o + 0.01833*Na_o)/(0.01833*Na_i +
        states[8])));
      const double cse_92 = PCa*cse_60*cse_74*states[0]*states[27]*states[28];
      const double cse_93 = Ca_o*std::exp(eta*cse_22)/(cse_2*cse_2*cse_2) +
        cse_1*cse_26*cse_59;
      const double cse_94 = K_SR*cse_85*cse_86;
      const double cse_95 =
        -0.5*G_KrMax*std::sqrt(K_o)*cse_58*states[4]/(1.4945*std::exp(0.0446*states[0])
        + 1.0);
      const double cse_96 = G_tiMax*cse_49*cse_58/((K_mK1 +
        K_o)*(std::pow(cse_0, -1.5)*std::exp(0.0561625551925629*states[0]) +
        2.0));
      const double cse_97 = k_NaCa*cse_41*cse_65*cse_84*cse_93;
      const double cse_98 =
        -3613.12438405488*PK*states[0]*states[27]*states[28]*(cse_34*states[8]
        + cse_49)/((1.0 + fmin(14452.4975362195*PCa*cse_60*cse_74*states[0],
        0)/ICahalf)*(cse_34 - 1.0));
      const double cse_99 = cse_87 + cse_90 + cse_91 + cse_95 + cse_96 +
        cse_98;
      values[0] =
        G_NaMax*cse_29*cse_61*states[2]*(states[3]*states[3]*states[3]) -
        G_bNaMax*cse_61 - ist + cse_26*cse_53 + cse_49*cse_89 - cse_79 -
        14452.4975362195*cse_92 - 5000.0*cse_97 + cse_99;
      values[1] = cse_29*(cse_83*(3.56*std::exp(0.079*states[0]) +
        310000.0*std::exp(0.35*states[0])) + (-1.0*cse_83 +
        1.0)/(0.13*std::exp(-0.0900900900900901*states[0] - 0.96036036036036)
        + 0.13)) + 0.135*cse_83*(cse_29 +
        1.0)*std::exp(-0.147058823529412*states[0] - 11.7647058823529);
      values[2] =
        cse_31*(0.1212*cse_83*std::exp(-0.01052*states[0])/(std::exp(-0.1378*states[0]
        - 5.531292) + 1.0) + (-0.3*cse_83 +
        0.3)*std::exp(-2.535e-7*states[0])/(std::exp(cse_36 - 3.2) + 1.0)) +
        cse_83*(cse_31 + 1.0)*(states[0] +
        37.78)*(-127140.0*std::exp(0.2444*states[0]) -
        3.474e-5*std::exp(-0.04391*states[0]))/(std::exp(0.311*states[0] +
        24.64053) + 1.0);
      values[3] = (1.0 - 1.0/(std::exp(cse_35 - 90.0) +
        1.0))*(-0.08*states[3]*std::exp(-0.0909090909090909*states[0]) +
        (-states[3] + 1.0)*(std::fabs(states[0] + 47.13) <= 1.0e-6 ?
        1.0/(-0.005*states[0] - 0.13565) : (0.32*states[0] +
        15.0816)/(-std::exp(cse_36 - 4.713) + 1.0)));
      values[4] = (cse_54*cse_88 - states[4])/(1.0*cse_88 + 27.0);
      values[5] = (-states[5] + 1.0/(std::exp(-0.0735294117647059*states[0] +
        1.81617647058824) + 1.0))*(1.0*(7.19e-5*states[0] -
        0.000719)/(-std::exp(-0.148*states[0] + 1.48) + 1.0) +
        1.0*(0.000131*states[0] - 0.00131)/(std::exp(0.0687*states[0] -
        0.687) - 1.0));
      values[6] = -0.0989*states[6]*std::exp(-0.06237*states[0]) +
        0.04516*(-states[6] + 1.0)*std::exp(0.03577*states[0]);
      values[7] = -0.005415*cse_55*states[7]/(0.051335*cse_55 + 1.0) +
        0.005415*cse_62*(cse_23 + 1.0)/(0.051335*cse_62 + 1.0);
      values[8] = 1.03626943005181e-5*cse_19*cse_6*(2.0*K_o*cse_89 + cse_99);
      values[9] = 1.0*(cse_50 - cse_73)/(CSQNtot*KmCSQN*std::pow(KmCSQN +
        states[9], -2.0) + 1.0);
      values[10] = -V_JSR*cse_3*cse_50 + V_myo*cse_3*cse_94;
      values[11] = 1.0*(-HTRPNtot*cse_76 - LTRPNtot*cse_77 -
        5.18134715025907e-6*cse_19*cse_6*(cse_53*states[11] + cse_79 -
        10000.0*cse_97) + cse_48 - cse_94)/(cse_21*std::pow(KmCMDN +
        states[11], -2.0) + 1.0);
      values[12] = 1.0*(V_JSR*cse_73*cse_9 - V_myo*cse_48*cse_9 -
        0.0748834069234172*cse_19*cse_9*cse_92)/(cse_21*std::pow(KmCMDN +
        states[12], -2.0) + 1.0);
      values[13] = kaminus*states[15] - cse_72;
      values[14] = kcplus*states[15] - cse_32;
      values[15] = kaminus*cse_24 + kcplus*cse_24 + cse_20 + cse_24*cse_57 +
        cse_32 + cse_72;
      values[16] = -cse_20 + cse_57*states[15];
      values[17] = omega*states[22] + cse_80*states[18] - states[17]*(cse_18
        + cse_67);
      values[18] = omega*cse_28 + cse_67*states[17] + cse_78*states[19] -
        states[18]*(aL*cse_18 + cse_71 + cse_80);
      values[19] = cse_51*states[24] + cse_71*states[18] + cse_82*states[20]
        - states[19]*(cse_14*cse_18 + cse_70 + cse_78);
      values[20] = cse_47*states[25] + cse_70*states[19] + cse_81*states[21]
        - states[20]*(cse_15*cse_18 + cse_69 + cse_82);
      values[21] = cse_27 + cse_52*states[26] + cse_69*states[20] -
        states[21]*(fL + cse_16*cse_18 + cse_81);
      values[22] = cse_18*states[17] + cse_28*cse_80 - states[22]*(aL*cse_67
        + omega);
      values[23] = aL*cse_18*states[18] + aL*cse_67*states[22] +
        cse_13*cse_78*states[24] - states[23]*(omega*cse_13 + cse_13*cse_80 +
        1.2*cse_68);
      values[24] = cse_13*cse_82*states[25] + cse_14*cse_18*states[19] +
        1.2*cse_68*states[23] - states[24]*(cse_13*cse_78 + cse_51 +
        0.8*cse_68);
      values[25] = cse_13*cse_81*states[26] + cse_15*cse_18*states[20] +
        0.8*cse_68*states[24] - states[25]*(cse_13*cse_82 + cse_47 +
        0.4*cse_68);
      values[26] = cse_16*cse_18*states[21] + 0.4*cse_68*states[25] -
        states[26]*(cse_13*cse_81 + cse_52);
      values[27] = fL*states[21] - cse_27;
      values[28] = (-states[28] + 0.2 + 0.8/(std::exp(cse_17 + 2.5) +
        1.0))/(20.0 + 600.0/(std::exp(0.105263157894737*states[0] +
        2.10526315789474) + 1.0));
      values[29] = cse_76;
      values[30] = cse_77;

    }

    // Get default initial conditions
    void get_ic(goss::DoubleVector *values) const
    {

      // Initial conditions
      values->n = _num_states;
      values->data.reset(new double[_num_states]);
      values->data[0] = -35;
      values->data[1] = 0.99869;
      values->data[2] = 0.99887;
      values->data[3] = 0.00024676;
      values->data[4] = 0.6935;
      values->data[5] = 0.00014589;
      values->data[6] = 3.742e-05;
      values->data[7] = 1;
      values->data[8] = 159.48;
      values->data[9] = 0.2616;
      values->data[10] = 0.262;
      values->data[11] = 8.464e-05;
      values->data[12] = 0.0001315;
      values->data[13] = 0.4929;
      values->data[14] = 0.5065;
      values->data[15] = 0.0006027;
      values->data[16] = 2.882e-09;
      values->data[17] = 0.99802;
      values->data[18] = 1.9544e-06;
      values->data[19] = 0;
      values->data[20] = 0;
      values->data[21] = 0;
      values->data[22] = 0.0019734;
      values->data[23] = 0;
      values->data[24] = 0;
      values->data[25] = 0;
      values->data[26] = 0;
      values->data[27] = 0;
      values->data[28] = 0.7959;
      values->data[29] = 0.13664;
      values->data[30] = 0.0055443;
    }

    // Return a copy of the ODE
    std::shared_ptr<ODE> copy() const
    {
      return std::make_shared<WinslowCSEArray>(*this);
    }

    // Evaluate the monitored intermediates
    void eval_monitored(const double* states, double t, double* monitored) const
    {

      // Assign states

      // Common Sub Expressions for monitored intermediates
      const double cse_0 = (Na_o*Na_o*Na_o);
      const double cse_1 = 0.0374417034617086*states[0];
      const double cse_2 = states[27]*states[28];
      const double cse_3 = std::exp(0.0748834069234172*states[0]);
      const double cse_4 = std::exp(cse_1);
      const double cse_5 = 1.0/(cse_3 - 1.0);
      const double cse_6 = std::exp(cse_1*(eta - 1.0));
      const double cse_7 = -0.341*Ca_o + 0.001*cse_3;
      const double cse_8 = 14452.4975362195*PCa*cse_5*cse_7*states[0];

      // Monitored intermediates
      monitored[1] = cse_2*cse_8;
      monitored[2] = 3613.12438405488*PK*cse_2*states[0]*(-K_o +
        cse_4*states[8])/((1.0 + fmin(cse_8, 0)/ICahalf)*(cse_4 - 1.0));
      monitored[3] = 5000.0*k_NaCa*(Ca_o*(Na_i*Na_i*Na_i)*std::exp(eta*cse_1)
        - cse_0*cse_6*states[11])/((Ca_o + K_mCa)*((K_mNa*K_mNa*K_mNa) +
        cse_0)*(k_sat*cse_6 + 1.0));
      monitored[4] = I_pCaMax*states[11]/(K_mpCa + states[11]);

    }

    // Set all field parameters
    void set_field_parameters(const double* field_params)
    {

      // Set field parameters
      G_KrMax = field_params[0];
      G_KsMax = field_params[1];

    }

  private:

    // Parameters
    double A_cap, C_sc, V_JSR, V_NSR, V_myo, V_ss, ist, Ca_o, K_o, Na_i,
      Na_o, G_KpMax, G_KrMax, G_KsMax, G_NaMax, G_bCaMax, G_bNaMax, G_tiMax,
      G_toMax, I_NaKMax, I_pCaMax, K_mCa, K_mK1, K_mKo, K_mNa, K_mNai,
      K_mpCa, eta, k_NaCa, k_sat, K_SR, K_fb, K_rb, N_fb, N_rb, kaminus,
      kaplus, kbminus, kbplus, kcminus, kcplus, mcoop, ncoop, tau_tr,
      tau_xfer, v_1, v_maxf, v_maxr, ICahalf, PCa, PK, aL, bL, fL, gL, omega,
      CMDNtot, CSQNtot, EGTAtot, HTRPNtot, KmCMDN, KmCSQN, KmEGTA, LTRPNtot,
      khtrpn_minus, khtrpn_plus, kltrpn_minus, kltrpn_plus;

  };

}

#endif
