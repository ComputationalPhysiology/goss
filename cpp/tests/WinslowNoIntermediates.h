#ifndef WINSLOWNOINTERMEDIATES_H_IS_INCLUDED
#define WINSLOWNOINTERMEDIATES_H_IS_INCLUDED

#include <stdexcept>
#include <cmath>
#include <memory>

#include "goss/ParameterizedODE.h"

namespace goss {

  // Implementation of gotran generated ODE
  class WinslowNoIntermediates : public ParameterizedODE
  {
  public:

    // Constructor
    WinslowNoIntermediates() : ParameterizedODE(31, 68, 1, 2, 4),
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
      const double V = states[0];
      const double h = states[1];
      const double j = states[2];
      const double m = states[3];
      const double xKr = states[4];
      const double xKs = states[5];
      const double xto1 = states[6];
      const double yto1 = states[7];
      const double K_i = states[8];
      const double Ca_JSR = states[9];
      const double Ca_NSR = states[10];
      const double Ca_i = states[11];
      const double Ca_ss = states[12];
      const double C1_RyR = states[13];
      const double C2_RyR = states[14];
      const double O1_RyR = states[15];
      const double O2_RyR = states[16];
      const double C0 = states[17];
      const double C1 = states[18];
      const double C2 = states[19];
      const double C3 = states[20];
      const double C4 = states[21];
      const double CCa0 = states[22];
      const double CCa1 = states[23];
      const double CCa2 = states[24];
      const double CCa3 = states[25];
      const double CCa4 = states[26];
      const double Open = states[27];
      const double yCa = states[28];
      const double HTRPNCa = states[29];
      const double LTRPNCa = states[30];
      values[0] = -Ca_i*I_pCaMax/(Ca_i + K_mpCa) - 1.0*G_KpMax*(V -
        26.7081865284974*std::log(K_o/K_i))/(std::exp(-0.167224080267559*V +
        1.25217391304348) + 1.0) - 0.5*G_KrMax*std::sqrt(K_o)*xKr*(V -
        26.7081865284974*std::log(K_o/K_i))/(1.4945*std::exp(0.0446*V) + 1.0)
        - G_KsMax*(xKs*xKs)*(V - 26.7081865284974*std::log((K_o +
        0.01833*Na_o)/(K_i + 0.01833*Na_i))) - G_NaMax*h*j*(m*m*m)*(V -
        26.7081865284974*std::log(Na_o/Na_i)) - G_bCaMax*(V -
        13.3540932642487*std::log(Ca_o/Ca_i)) - G_bNaMax*(V -
        26.7081865284974*std::log(Na_o/Na_i)) - 1.0*G_tiMax*K_o*(V -
        26.7081865284974*std::log(K_o/K_i))/((K_mK1 + K_o)*(std::pow(K_o/K_i,
        -1.5)*std::exp(0.0561625551925629*V) + 2.0)) - G_toMax*xto1*yto1*(V -
        26.7081865284974*std::log(K_o/K_i)) - 1.0*I_NaKMax*K_o/((1.0 +
        0.1245*std::exp(-0.00374417034617086*V))*(K_mKo +
        K_o)*(std::pow(K_mNai/Na_i, 1.5) + 1.0)) -
        14452.4975362195*Open*PCa*V*yCa*(-0.341*Ca_o +
        0.001*std::exp(0.0748834069234172*V))/(std::exp(0.0748834069234172*V)
        - 1.0) -
        3613.12438405488*Open*PK*V*yCa*(K_i*std::exp(0.0374417034617086*V) -
        K_o)/((1.0 + fmin(14452.4975362195*PCa*V*(-0.341*Ca_o +
        0.001*std::exp(0.0748834069234172*V))/(std::exp(0.0748834069234172*V)
        - 1.0), 0)/ICahalf)*(std::exp(0.0374417034617086*V) - 1.0)) - ist -
        5000.0*k_NaCa*(-Ca_i*(Na_o*Na_o*Na_o)*std::exp(0.0374417034617086*V*(eta
        - 1.0)) +
        Ca_o*(Na_i*Na_i*Na_i)*std::exp(0.0374417034617086*V*eta))/((Ca_o +
        K_mCa)*((K_mNa*K_mNa*K_mNa) +
        (Na_o*Na_o*Na_o))*(k_sat*std::exp(0.0374417034617086*V*(eta - 1.0)) +
        1.0));
      values[1] = -h*((1.0 - 1.0/(std::exp(-1.0*V - 40.0) +
        1.0))/(0.13*std::exp(-0.0900900900900901*V - 0.96036036036036) +
        0.13) + (3.56*std::exp(0.079*V) +
        310000.0*std::exp(0.35*V))/(std::exp(-1.0*V - 40.0) + 1.0)) +
        0.135*(-h + 1.0)*std::exp(-0.147058823529412*V -
        11.7647058823529)/(std::exp(-1.0*V - 40.0) + 1.0);
      values[2] = -j*((0.3 - 0.3/(std::exp(-1.0*V - 40.0) +
        1.0))*std::exp(-2.535e-7*V)/(std::exp(-0.1*V - 3.2) + 1.0) +
        0.1212*std::exp(-0.01052*V)/((std::exp(-1.0*V - 40.0) +
        1.0)*(std::exp(-0.1378*V - 5.531292) + 1.0))) + (V + 37.78)*(-j +
        1.0)*(-127140.0*std::exp(0.2444*V) -
        3.474e-5*std::exp(-0.04391*V))/((std::exp(-1.0*V - 40.0) +
        1.0)*(std::exp(0.311*V + 24.64053) + 1.0));
      values[3] = (1.0 - 1.0/(std::exp(-1.0*V - 90.0) +
        1.0))*(-0.08*m*std::exp(-0.0909090909090909*V) + (-m +
        1.0)*(std::fabs(V + 47.13) <= 1.0e-6 ? 1.0/(-0.005*V - 0.13565) :
        (0.32*V + 15.0816)/(-std::exp(-0.1*V - 4.713) + 1.0)));
      values[4] = (-xKr + std::exp(0.1691*V - 5.495)/(std::exp(-0.0128*V -
        7.677) + std::exp(0.1691*V - 5.495)))/(27.0 + 1.0/(std::exp(-0.0128*V
        - 7.677) + std::exp(0.1691*V - 5.495)));
      values[5] = (-xKs + 1.0/(std::exp(-0.0735294117647059*V +
        1.81617647058824) + 1.0))*(1.0*(7.19e-5*V -
        0.000719)/(-std::exp(-0.148*V + 1.48) + 1.0) + 1.0*(0.000131*V -
        0.00131)/(std::exp(0.0687*V - 0.687) - 1.0));
      values[6] = -0.0989*xto1*std::exp(-0.06237*V) + 0.04516*(-xto1 +
        1.0)*std::exp(0.03577*V);
      values[7] = -0.005415*yto1*std::exp(0.2*V +
        6.7)/(0.051335*std::exp(0.2*V + 6.7) + 1.0) + 0.005415*(-yto1 +
        1.0)*std::exp(-0.2*V - 6.7)/(0.051335*std::exp(-0.2*V - 6.7) + 1.0);
      values[8] = 1.03626943005181e-5*A_cap*C_sc*(-1.0*G_KpMax*(V -
        26.7081865284974*std::log(K_o/K_i))/(std::exp(-0.167224080267559*V +
        1.25217391304348) + 1.0) - 0.5*G_KrMax*std::sqrt(K_o)*xKr*(V -
        26.7081865284974*std::log(K_o/K_i))/(1.4945*std::exp(0.0446*V) + 1.0)
        - G_KsMax*(xKs*xKs)*(V - 26.7081865284974*std::log((K_o +
        0.01833*Na_o)/(K_i + 0.01833*Na_i))) - 1.0*G_tiMax*K_o*(V -
        26.7081865284974*std::log(K_o/K_i))/((K_mK1 + K_o)*(std::pow(K_o/K_i,
        -1.5)*std::exp(0.0561625551925629*V) + 2.0)) - G_toMax*xto1*yto1*(V -
        26.7081865284974*std::log(K_o/K_i)) + 2.0*I_NaKMax*K_o/((1.0 +
        0.1245*std::exp(-0.00374417034617086*V))*(K_mKo +
        K_o)*(std::pow(K_mNai/Na_i, 1.5) + 1.0)) -
        3613.12438405488*Open*PK*V*yCa*(K_i*std::exp(0.0374417034617086*V) -
        K_o)/((1.0 + fmin(14452.4975362195*PCa*V*(-0.341*Ca_o +
        0.001*std::exp(0.0748834069234172*V))/(std::exp(0.0748834069234172*V)
        - 1.0), 0)/ICahalf)*(std::exp(0.0374417034617086*V) - 1.0)))/V_myo;
      values[9] = 1.0*(-v_1*(Ca_JSR - Ca_ss)*(O1_RyR + O2_RyR) + (-Ca_JSR +
        Ca_NSR)/tau_tr)/(CSQNtot*KmCSQN*std::pow(Ca_JSR + KmCSQN, -2.0) +
        1.0);
      values[10] = K_SR*V_myo*(v_maxf*std::pow(Ca_i/K_fb, N_fb) -
        v_maxr*std::pow(Ca_NSR/K_rb, N_rb))/(V_NSR*(std::pow(Ca_NSR/K_rb,
        N_rb) + std::pow(Ca_i/K_fb, N_fb) + 1.0)) - V_JSR*(-Ca_JSR +
        Ca_NSR)/(V_NSR*tau_tr);
      values[11] = 1.0*(-5.18134715025907e-6*A_cap*C_sc*(Ca_i*I_pCaMax/(Ca_i
        + K_mpCa) + G_bCaMax*(V - 13.3540932642487*std::log(Ca_o/Ca_i)) -
        10000.0*k_NaCa*(-Ca_i*(Na_o*Na_o*Na_o)*std::exp(0.0374417034617086*V*(eta
        - 1.0)) +
        Ca_o*(Na_i*Na_i*Na_i)*std::exp(0.0374417034617086*V*eta))/((Ca_o +
        K_mCa)*((K_mNa*K_mNa*K_mNa) +
        (Na_o*Na_o*Na_o))*(k_sat*std::exp(0.0374417034617086*V*(eta - 1.0)) +
        1.0)))/V_myo - HTRPNtot*(Ca_i*khtrpn_plus*(-HTRPNCa + 1.0) -
        HTRPNCa*khtrpn_minus) - K_SR*(v_maxf*std::pow(Ca_i/K_fb, N_fb) -
        v_maxr*std::pow(Ca_NSR/K_rb, N_rb))/(std::pow(Ca_NSR/K_rb, N_rb) +
        std::pow(Ca_i/K_fb, N_fb) + 1.0) -
        LTRPNtot*(Ca_i*kltrpn_plus*(-LTRPNCa + 1.0) - LTRPNCa*kltrpn_minus) +
        (-Ca_i + Ca_ss)/tau_xfer)/(CMDNtot*KmCMDN*std::pow(Ca_i + KmCMDN,
        -2.0) + 1.0);
      values[12] =
        1.0*(-0.0748834069234172*A_cap*C_sc*Open*PCa*V*yCa*(-0.341*Ca_o +
        0.001*std::exp(0.0748834069234172*V))/(V_ss*(std::exp(0.0748834069234172*V)
        - 1.0)) + V_JSR*v_1*(Ca_JSR - Ca_ss)*(O1_RyR + O2_RyR)/V_ss -
        V_myo*(-Ca_i + Ca_ss)/(V_ss*tau_xfer))/(CMDNtot*KmCMDN*std::pow(Ca_ss
        + KmCMDN, -2.0) + 1.0);
      values[13] = -C1_RyR*kaplus*std::pow(1000.0*Ca_ss, ncoop) +
        O1_RyR*kaminus;
      values[14] = -C2_RyR*kcminus + O1_RyR*kcplus;
      values[15] = C1_RyR*kaplus*std::pow(1000.0*Ca_ss, ncoop) +
        C2_RyR*kcminus - O1_RyR*kaminus -
        O1_RyR*kbplus*std::pow(1000.0*Ca_ss, mcoop) - O1_RyR*kcplus +
        O2_RyR*kbminus;
      values[16] = O1_RyR*kbplus*std::pow(1000.0*Ca_ss, mcoop) -
        O2_RyR*kbminus;
      values[17] = -C0*(0.10375*Ca_ss + 1.6*std::exp(0.1*V + 0.2)) +
        0.05*C1*std::exp(-0.0769230769230769*V - 0.153846153846154) +
        CCa0*omega;
      values[18] = 1.6*C0*std::exp(0.1*V + 0.2) - C1*(0.10375*Ca_ss*aL +
        0.05*std::exp(-0.0769230769230769*V - 0.153846153846154) +
        1.2*std::exp(0.1*V + 0.2)) + 0.1*C2*std::exp(-0.0769230769230769*V -
        0.153846153846154) + CCa1*omega/bL;
      values[19] = 1.2*C1*std::exp(0.1*V + 0.2) - C2*(0.10375*Ca_ss*(aL*aL) +
        0.1*std::exp(-0.0769230769230769*V - 0.153846153846154) +
        0.8*std::exp(0.1*V + 0.2)) + 0.15*C3*std::exp(-0.0769230769230769*V -
        0.153846153846154) + CCa2*omega/(bL*bL);
      values[20] = 0.8*C2*std::exp(0.1*V + 0.2) -
        C3*(0.10375*Ca_ss*(aL*aL*aL) + 0.15*std::exp(-0.0769230769230769*V -
        0.153846153846154) + 0.4*std::exp(0.1*V + 0.2)) +
        0.2*C4*std::exp(-0.0769230769230769*V - 0.153846153846154) +
        CCa3*omega/(bL*bL*bL);
      values[21] = 0.4*C3*std::exp(0.1*V + 0.2) -
        C4*(0.10375*Ca_ss*(aL*aL*aL*aL) + fL +
        0.2*std::exp(-0.0769230769230769*V - 0.153846153846154)) +
        CCa4*omega/(bL*bL*bL*bL) + Open*gL;
      values[22] = 0.10375*C0*Ca_ss - CCa0*(1.6*aL*std::exp(0.1*V + 0.2) +
        omega) + 0.05*CCa1*std::exp(-0.0769230769230769*V -
        0.153846153846154)/bL;
      values[23] = 0.10375*C1*Ca_ss*aL + 1.6*CCa0*aL*std::exp(0.1*V + 0.2) -
        CCa1*(1.2*aL*std::exp(0.1*V + 0.2) + omega/bL +
        0.05*std::exp(-0.0769230769230769*V - 0.153846153846154)/bL) +
        0.1*CCa2*std::exp(-0.0769230769230769*V - 0.153846153846154)/bL;
      values[24] = 0.10375*C2*Ca_ss*(aL*aL) + 1.2*CCa1*aL*std::exp(0.1*V +
        0.2) - CCa2*(0.8*aL*std::exp(0.1*V + 0.2) +
        0.1*std::exp(-0.0769230769230769*V - 0.153846153846154)/bL +
        omega/(bL*bL)) + 0.15*CCa3*std::exp(-0.0769230769230769*V -
        0.153846153846154)/bL;
      values[25] = 0.10375*C3*Ca_ss*(aL*aL*aL) + 0.8*CCa2*aL*std::exp(0.1*V +
        0.2) - CCa3*(0.4*aL*std::exp(0.1*V + 0.2) +
        0.15*std::exp(-0.0769230769230769*V - 0.153846153846154)/bL +
        omega/(bL*bL*bL)) + 0.2*CCa4*std::exp(-0.0769230769230769*V -
        0.153846153846154)/bL;
      values[26] = 0.10375*C4*Ca_ss*(aL*aL*aL*aL) +
        0.4*CCa3*aL*std::exp(0.1*V + 0.2) -
        CCa4*(0.2*std::exp(-0.0769230769230769*V - 0.153846153846154)/bL +
        omega/(bL*bL*bL*bL));
      values[27] = C4*fL - Open*gL;
      values[28] = (-yCa + 0.2 + 0.8/(std::exp(0.2*V + 2.5) + 1.0))/(20.0 +
        600.0/(std::exp(0.105263157894737*V + 2.10526315789474) + 1.0));
      values[29] = Ca_i*khtrpn_plus*(-HTRPNCa + 1.0) - HTRPNCa*khtrpn_minus;
      values[30] = Ca_i*kltrpn_plus*(-LTRPNCa + 1.0) - LTRPNCa*kltrpn_minus;

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
      return std::make_shared<WinslowNoIntermediates>(*this);
    }

    // Evaluate the monitored intermediates
    void eval_monitored(const double* states, double t, double* monitored) const
    {

      // Assign states
      const double V = states[0];
      const double K_i = states[8];
      const double Ca_i = states[11];
      const double Open = states[27];
      const double yCa = states[28];

      // Common Sub Expressions for monitored intermediates
      const double cse_0 = (Na_o*Na_o*Na_o);
      const double cse_1 = 0.0374417034617086*V;
      const double cse_2 = std::exp(0.0748834069234172*V);
      const double cse_3 = std::exp(cse_1);
      const double cse_4 = Open*V*yCa;
      const double cse_5 = 1.0/(cse_2 - 1.0);
      const double cse_6 = std::exp(cse_1*(eta - 1.0));
      const double cse_7 = -0.341*Ca_o + 0.001*cse_2;
      const double cse_8 = 14452.4975362195*PCa*cse_5*cse_7;

      // Monitored intermediates
      monitored[1] = cse_4*cse_8;
      monitored[2] = 3613.12438405488*PK*cse_4*(K_i*cse_3 - K_o)/((1.0 +
        fmin(V*cse_8, 0)/ICahalf)*(cse_3 - 1.0));
      monitored[3] = 5000.0*k_NaCa*(-Ca_i*cse_0*cse_6 +
        Ca_o*(Na_i*Na_i*Na_i)*std::exp(eta*cse_1))/((Ca_o +
        K_mCa)*((K_mNa*K_mNa*K_mNa) + cse_0)*(k_sat*cse_6 + 1.0));
      monitored[4] = Ca_i*I_pCaMax/(Ca_i + K_mpCa);

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
