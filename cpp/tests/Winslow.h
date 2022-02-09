#ifndef WINSLOW_H_IS_INCLUDED
#define WINSLOW_H_IS_INCLUDED

#include <stdexcept>
#include <cmath>
#include <memory>

#include "goss/ParameterizedODE.h"

namespace goss {

  // Implementation of gotran generated ODE
  class Winslow : public ParameterizedODE
  {
  public:

    // Constructor
    Winslow() : ParameterizedODE(31, 68, 1, 2, 4),
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

      // Constants
      const double F = 96.5;
      const double T = 310;
      const double R = 8.314;
      const double RTonF = R*T/F;
      const double FonRT = F/(R*T);

      // Help variables
      const double VFonRT = FonRT*V;
      const double expVFonRT = std::exp(VFonRT);

      // I Membrane currents

      // Na+ current I_Na
      const double E_Na = RTonF*std::log(Na_o/Na_i);
      const double I_Na = G_NaMax*h*j*(m*m*m)*(-E_Na + V);
      const double a_h = 0.135*std::exp(-0.147058823529412*V -
        11.7647058823529)/(std::exp(-1.0*V - 40.0) + 1.0);
      const double b_h = (1.0 - 1.0/(std::exp(-1.0*V - 40.0) +
        1.0))/(0.13*std::exp(-0.0900900900900901*V - 0.96036036036036) +
        0.13) + (3.56*std::exp(0.079*V) +
        310000.0*std::exp(0.35*V))/(std::exp(-1.0*V - 40.0) + 1.0);
      const double a_j = (V + 37.78)*(-127140.0*std::exp(0.2444*V) -
        3.474e-5*std::exp(-0.04391*V))/((std::exp(-1.0*V - 40.0) +
        1.0)*(std::exp(0.311*V + 24.64053) + 1.0));
      const double b_j = 0.3*(1.0 - 1.0/(std::exp(-1.0*V - 40.0) +
        1.0))*std::exp(-2.535e-7*V)/(std::exp(-0.1*V - 3.2) + 1.0) +
        0.1212*std::exp(-0.01052*V)/((std::exp(-1.0*V - 40.0) +
        1.0)*(std::exp(-0.1378*V - 5.531292) + 1.0));
      const double a_m = (std::fabs(V + 47.13) <= 1.0e-6 ? 1.0/(-0.005*V -
        0.13565) : (0.32*V + 15.0816)/(-std::exp(-0.1*V - 4.713) + 1.0));
      const double b_m = 0.08*std::exp(-V/11.0);
      const double dm = (1.0 - 1.0/(std::exp(-1.0*V - 90.0) + 1.0))*(a_m*(-m
        + 1.0) - b_m*m);

      // Rapid-activating delayed rectifier K+ current I_Kr
      const double k12 = std::exp(0.1691*V - 5.495);
      const double k21 = std::exp(-0.0128*V - 7.677);
      const double xKr_inf = k12/(k12 + k21);
      const double tau_xKr = 27.0 + 1.0/(k12 + k21);
      const double dxKr = (-xKr + xKr_inf)/tau_xKr;
      const double E_k = RTonF*std::log(K_o/K_i);
      const double R_V = 1.0/(1.4945*std::exp(0.0446*V) + 1.0);
      const double f_k = std::sqrt(K_o)/2.0;
      const double I_Kr = G_KrMax*R_V*f_k*xKr*(-E_k + V);

      // Slow-activating delayed rectifier K+ current I_Ks
      const double xKs_inf = 1.0/(std::exp(-0.0735294117647059*V +
        1.81617647058824) + 1.0);
      const double tau_xKs = 1.0/((7.19e-5*V - 0.000719)/(-std::exp(-0.148*V
        + 1.48) + 1.0) + (0.000131*V - 0.00131)/(std::exp(0.0687*V - 0.687) -
        1.0));
      const double dxKs = (-xKs + xKs_inf)/tau_xKs;
      const double E_Ks = RTonF*std::log((K_o + 0.01833*Na_o)/(K_i +
        0.01833*Na_i));
      const double I_Ks = G_KsMax*(xKs*xKs)*(-E_Ks + V);

      // Transient outward K+ current I_to
      const double alpha_xto1 = 0.04516*std::exp(0.03577*V);
      const double beta_xto1 = 0.0989*std::exp(-0.06237*V);
      double a1 = 0.051335*std::exp(-0.2*V - 6.7) + 1.0;
      const double alpha_yto1 = 0.005415*std::exp(-0.2*V - 6.7)/a1;
      a1 = 0.051335*std::exp(0.2*V + 6.7) + 1.0;
      const double beta_yto1 = 0.005415*std::exp(0.2*V + 6.7)/a1;
      const double dxto1 = alpha_xto1*(-xto1 + 1.0) - beta_xto1*xto1;
      const double dyto1 = alpha_yto1*(-yto1 + 1.0) - beta_yto1*yto1;
      const double I_to = G_toMax*xto1*yto1*(-E_k + V);

      // Time-Independent K+ current I_ti
      const double K_tiUnlim = 1.0/(std::exp(FonRT*(-1.5*E_k + 1.5*V)) + 2.0);
      const double I_ti = G_tiMax*K_o*K_tiUnlim*(-E_k + V)/(K_mK1 + K_o);

      // Plateau current I_Kp
      const double K_p = 1.0/(std::exp(-0.167224080267559*V +
        1.25217391304348) + 1.0);
      const double I_Kp = G_KpMax*K_p*(-E_k + V);

      // NCX Current I_NaCa
      const double I_NaCa =
        5000.0*k_NaCa*(-Ca_i*(Na_o*Na_o*Na_o)*std::exp(VFonRT*(eta - 1.0)) +
        Ca_o*(Na_i*Na_i*Na_i)*std::exp(VFonRT*eta))/((Ca_o +
        K_mCa)*((K_mNa*K_mNa*K_mNa) +
        (Na_o*Na_o*Na_o))*(k_sat*std::exp(FonRT*V*(eta - 1.0)) + 1.0));

      // Na+-K+ pump current I_NaK
      const double sigma = 0;
      const double f_NaK = 1.0/(0.0365*sigma*std::exp(-VFonRT) + 1.0 +
        0.1245*std::exp(-0.1*VFonRT));
      const double I_NaK = I_NaKMax*K_o*f_NaK/((K_mKo +
        K_o)*(std::pow(K_mNai/Na_i, 1.5) + 1.0));

      // Sarcolemmal Ca2+ pump current I_pCa
      const double I_pCa = Ca_i*I_pCaMax/(Ca_i + K_mpCa);

      // Ca2+ background current I_bCa
      const double E_Ca = RTonF*std::log(Ca_o/Ca_i)/2.0;
      const double I_bCa = G_bCaMax*(-E_Ca + V);

      // Na+ background current I_bNa
      const double I_bNa = G_bNaMax*(-E_Na + V);

      // II Ca2+ handling mechanisms

      // L-type Ca2+ current I_Ca
      const double alpha = 0.4*std::exp(0.1*V + 0.2);
      const double beta = 0.05*std::exp(-0.0769230769230769*V -
        0.153846153846154);
      const double alpha_prime = aL*alpha;
      const double beta_prime = beta/bL;
      double gamma = 0.10375*Ca_ss;
      const double C0_to_C1 = 4.0*alpha;
      const double C1_to_C2 = 3.0*alpha;
      const double C2_to_C3 = 2.0*alpha;
      const double C3_to_C4 = alpha;
      const double CCa0_to_CCa1 = 4.0*alpha_prime;
      const double CCa1_to_CCa2 = 3.0*alpha_prime;
      const double CCa2_to_CCa3 = 2.0*alpha_prime;
      const double CCa3_to_CCa4 = alpha_prime;
      const double C1_to_C0 = beta;
      const double C2_to_C1 = 2.0*beta;
      const double C3_to_C2 = 3.0*beta;
      const double C4_to_C3 = 4.0*beta;
      const double CCa1_to_CCa0 = beta_prime;
      const double CCa2_to_CCa1 = 2.0*beta_prime;
      const double CCa3_to_CCa2 = 3.0*beta_prime;
      const double CCa4_to_CCa3 = 4.0*beta_prime;
      gamma = 0.10375*Ca_ss;
      const double C0_to_CCa0 = gamma;
      const double C1_to_CCa1 = C0_to_CCa0*aL;
      const double C2_to_CCa2 = C1_to_CCa1*aL;
      const double C3_to_CCa3 = C2_to_CCa2*aL;
      const double C4_to_CCa4 = C3_to_CCa3*aL;
      const double CCa0_to_C0 = omega;
      const double CCa1_to_C1 = CCa0_to_C0/bL;
      const double CCa2_to_C2 = CCa1_to_C1/bL;
      const double CCa3_to_C3 = CCa2_to_C2/bL;
      const double CCa4_to_C4 = CCa3_to_C3/bL;
      a1 = C0*(C0_to_C1 + C0_to_CCa0);
      double a2 = C1*C1_to_C0 + CCa0*CCa0_to_C0;
      const double dC0 = -a1 + a2;
      a1 = C1*(C1_to_C0 + C1_to_C2 + C1_to_CCa1);
      a2 = C0*C0_to_C1 + C2*C2_to_C1 + CCa1*CCa1_to_C1;
      const double dC1 = -a1 + a2;
      a1 = C2*(C2_to_C1 + C2_to_C3 + C2_to_CCa2);
      a2 = C1*C1_to_C2 + C3*C3_to_C2 + CCa2*CCa2_to_C2;
      const double dC2 = -a1 + a2;
      a1 = C3*(C3_to_C2 + C3_to_C4 + C3_to_CCa3);
      a2 = C2*C2_to_C3 + C4*C4_to_C3 + CCa3*CCa3_to_C3;
      const double dC3 = -a1 + a2;
      a1 = C4*(C4_to_C3 + C4_to_CCa4 + fL);
      a2 = C3*C3_to_C4 + CCa4*CCa4_to_C4 + Open*gL;
      const double dC4 = -a1 + a2;
      const double dOpen = C4*fL - Open*gL;
      a1 = CCa0*(CCa0_to_C0 + CCa0_to_CCa1);
      a2 = C0*C0_to_CCa0 + CCa1*CCa1_to_CCa0;
      const double dCCa0 = -a1 + a2;
      a1 = CCa1*(CCa1_to_C1 + CCa1_to_CCa0 + CCa1_to_CCa2);
      a2 = C1*C1_to_CCa1 + CCa0*CCa0_to_CCa1 + CCa2*CCa2_to_CCa1;
      const double dCCa1 = -a1 + a2;
      a1 = CCa2*(CCa2_to_C2 + CCa2_to_CCa1 + CCa2_to_CCa3);
      a2 = C2*C2_to_CCa2 + CCa1*CCa1_to_CCa2 + CCa3*CCa3_to_CCa2;
      const double dCCa2 = -a1 + a2;
      a1 = CCa3*(CCa3_to_C3 + CCa3_to_CCa2 + CCa3_to_CCa4);
      a2 = C3*C3_to_CCa3 + CCa2*CCa2_to_CCa3 + CCa4*CCa4_to_CCa3;
      const double dCCa3 = -a1 + a2;
      a1 = CCa4*(CCa4_to_C4 + CCa4_to_CCa3);
      a2 = C4*C4_to_CCa4 + CCa3*CCa3_to_CCa4;
      const double dCCa4 = -a1 + a2;
      const double yCa_inf = 0.2 + 0.8/(std::exp(0.2*V + 2.5) + 1.0);
      const double tau_yCa = 20.0 + 600.0/(std::exp(0.105263157894737*V +
        2.10526315789474) + 1.0);
      const double dyCa = (-yCa + yCa_inf)/tau_yCa;
      const double VFsqonRT = 1000.0*F*VFonRT;
      a1 = -0.341*Ca_o + 0.001*std::exp(2.0*VFonRT);
      a2 = std::exp(2.0*VFonRT) - 1.0;
      const double ICamax = 4.0*PCa*VFsqonRT*a1/a2;
      const double I_Ca = ICamax*Open*yCa;
      const double PKprime = PK/(1.0 + fmin(ICamax, 0)/ICahalf);
      a1 = K_i*expVFonRT - K_o;
      a2 = expVFonRT - 1.0;
      const double I_CaK = Open*PKprime*VFsqonRT*a1*yCa/a2;

      // RyR Channel
      a1 = std::pow(1000.0*Ca_ss, mcoop);
      a2 = std::pow(1000.0*Ca_ss, ncoop);
      const double dC1_RyR = -C1_RyR*a2*kaplus + O1_RyR*kaminus;
      const double dO2_RyR = O1_RyR*a1*kbplus - O2_RyR*kbminus;
      const double dC2_RyR = -C2_RyR*kcminus + O1_RyR*kcplus;
      const double dO1_RyR = -dC1_RyR - dC2_RyR - dO2_RyR;
      const double J_rel = v_1*(Ca_JSR - Ca_ss)*(O1_RyR + O2_RyR);

      // SERCA2a Pump
      const double f_b = std::pow(Ca_i/K_fb, N_fb);
      const double r_b = std::pow(Ca_NSR/K_rb, N_rb);
      const double J_up = K_SR*(f_b*v_maxf - r_b*v_maxr)/(f_b + r_b + 1.0);

      // Intracellular Ca fluxes
      const double J_tr = (-Ca_JSR + Ca_NSR)/tau_tr;
      const double J_xfer = (-Ca_i + Ca_ss)/tau_xfer;
      a1 = LTRPNCa*kltrpn_minus;
      const double dLTRPNCa = Ca_i*kltrpn_plus*(-LTRPNCa + 1.0) - a1;
      a1 = HTRPNCa*khtrpn_minus;
      const double dHTRPNCa = Ca_i*khtrpn_plus*(-HTRPNCa + 1.0) - a1;
      const double J_trpn = HTRPNtot*dHTRPNCa + LTRPNtot*dLTRPNCa;
      a1 = CMDNtot*KmCMDN*std::pow(Ca_ss + KmCMDN, -2.0);
      const double beta_ss = 1.0/(a1 + 1.0);
      a1 = CSQNtot*KmCSQN*std::pow(Ca_JSR + KmCSQN, -2.0);
      const double beta_JSR = 1.0/(a1 + 1.0);
      a1 = CMDNtot*KmCMDN*std::pow(Ca_i + KmCMDN, -2.0);
      const double beta_i = 1.0/(a1 + 1.0);

      // The ODE system
      values[0] = -I_Ca - I_CaK - I_Kp - I_Kr - I_Ks - I_Na - I_NaCa - I_NaK
        - I_bCa - I_bNa - I_pCa - I_ti - I_to - ist;
      values[1] = a_h*(-h + 1.0) - b_h*h;
      values[2] = a_j*(-j + 1.0) - b_j*j;
      values[3] = dm;
      values[4] = dxKr;
      values[5] = dxKs;
      values[6] = dxto1;
      values[7] = dyto1;
      values[8] = A_cap*C_sc*(-I_CaK - I_Kp - I_Kr - I_Ks + 2.0*I_NaK - I_ti
        - I_to)/(1000.0*F*V_myo);
      values[9] = beta_JSR*(-J_rel + J_tr);
      values[10] = -J_tr*V_JSR/V_NSR + J_up*V_myo/V_NSR;
      values[11] = beta_i*(-A_cap*C_sc*(-2.0*I_NaCa + I_bCa +
        I_pCa)/(2000.0*F*V_myo) - J_trpn - J_up + J_xfer);
      values[12] = beta_ss*(-A_cap*C_sc*I_Ca/(2000.0*F*V_ss) +
        J_rel*V_JSR/V_ss - J_xfer*V_myo/V_ss);
      values[13] = dC1_RyR;
      values[14] = dC2_RyR;
      values[15] = dO1_RyR;
      values[16] = dO2_RyR;
      values[17] = dC0;
      values[18] = dC1;
      values[19] = dC2;
      values[20] = dC3;
      values[21] = dC4;
      values[22] = dCCa0;
      values[23] = dCCa1;
      values[24] = dCCa2;
      values[25] = dCCa3;
      values[26] = dCCa4;
      values[27] = dOpen;
      values[28] = dyCa;
      values[29] = dHTRPNCa;
      values[30] = dLTRPNCa;

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
      return std::make_shared<Winslow>(*this);
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
