import math

import benchmark
import numpy as np


def rhs(t, states, parameters):
    """
    Compute the right hand side of the tentusscher_panfilov_2006_M_cell ODE
    """

    # Assign states
    (
        Xr1,
        Xr2,
        Xs,
        m,
        h,
        j,
        d,
        f,
        f2,
        fCass,
        s,
        r,
        Ca_i,
        R_prime,
        Ca_SR,
        Ca_ss,
        Na_i,
        V,
        K_i,
    ) = states

    # Assign parameters
    (
        P_kna,
        g_K1,
        g_Kr,
        g_Ks,
        g_Na,
        g_bna,
        g_CaL,
        g_bca,
        g_to,
        K_mNa,
        K_mk,
        P_NaK,
        K_NaCa,
        K_sat,
        Km_Ca,
        Km_Nai,
        alpha,
        gamma,
        K_pCa,
        g_pCa,
        g_pK,
        Buf_c,
        Buf_sr,
        Buf_ss,
        Ca_o,
        EC,
        K_buf_c,
        K_buf_sr,
        K_buf_ss,
        K_up,
        V_leak,
        V_rel,
        V_sr,
        V_ss,
        V_xfer,
        Vmax_up,
        k1_prime,
        k2_prime,
        k3,
        k4,
        max_sr,
        min_sr,
        Na_o,
        Cm,
        F,
        R,
        T,
        V_c,
        stim_amplitude,
        stim_duration,
        stim_period,
        stim_start,
        K_o,
    ) = parameters

    # Init return args
    values = np.zeros((19,), dtype=np.float_)

    # Expressions for the Reversal potentials component
    E_Na = R * T * math.log(Na_o / Na_i) / F
    E_K = R * T * math.log(K_o / K_i) / F
    E_Ks = R * T * math.log((K_o + Na_o * P_kna) / (P_kna * Na_i + K_i)) / F
    E_Ca = 0.5 * R * T * math.log(Ca_o / Ca_i) / F

    # Expressions for the Inward rectifier potassium current component
    alpha_K1 = 0.1 / (1 + 6.14421235332821e-06 * math.exp(0.06 * V - 0.06 * E_K))
    beta_K1 = (
        0.36787944117144233 * math.exp(0.1 * V - 0.1 * E_K)
        + 3.0606040200802673 * math.exp(0.0002 * V - 0.0002 * E_K)
    ) / (1 + math.exp(0.5 * E_K - 0.5 * V))
    xK1_inf = alpha_K1 / (alpha_K1 + beta_K1)
    i_K1 = 0.4303314829119352 * g_K1 * math.sqrt(K_o) * (-E_K + V) * xK1_inf

    # Expressions for the Rapid time dependent potassium current component
    i_Kr = 0.4303314829119352 * g_Kr * math.sqrt(K_o) * (-E_K + V) * Xr1 * Xr2

    # Expressions for the Xr1 gate component
    xr1_inf = 1.0 / (1 + math.exp(-26 / 7 - V / 7))
    alpha_xr1 = 450 / (1 + math.exp(-9 / 2 - V / 10))
    beta_xr1 = 6 / (1 + 13.581324522578193 * math.exp(0.08695652173913043 * V))
    tau_xr1 = alpha_xr1 * beta_xr1
    values[0] = (-Xr1 + xr1_inf) / tau_xr1

    # Expressions for the Xr2 gate component
    xr2_inf = 1.0 / (1 + math.exp(11 / 3 + V / 24))
    alpha_xr2 = 3 / (1 + math.exp(-3 - V / 20))
    beta_xr2 = 1.12 / (1 + math.exp(-3 + V / 20))
    tau_xr2 = alpha_xr2 * beta_xr2
    values[1] = (-Xr2 + xr2_inf) / tau_xr2

    # Expressions for the Slow time dependent potassium current component
    i_Ks = g_Ks * (Xs * Xs) * (-E_Ks + V)

    # Expressions for the Xs gate component
    xs_inf = 1.0 / (1 + math.exp(-5 / 14 - V / 14))
    alpha_xs = 1400 / math.sqrt(1 + math.exp(5 / 6 - V / 6))
    beta_xs = 1.0 / (1 + math.exp(-7 / 3 + V / 15))
    tau_xs = 80 + alpha_xs * beta_xs
    values[2] = (-Xs + xs_inf) / tau_xs

    # Expressions for the Fast sodium current component
    i_Na = g_Na * (m * m * m) * (-E_Na + V) * h * j

    # Expressions for the m gate component
    m_inf = 1.0 / (
        (1 + 0.0018422115811651339 * math.exp(-0.1107419712070875 * V))
        * (1 + 0.0018422115811651339 * math.exp(-0.1107419712070875 * V))
    )
    alpha_m = 1.0 / (1 + math.exp(-12 - V / 5))
    beta_m = 0.1 / (1 + math.exp(7 + V / 5)) + 0.1 / (1 + math.exp(-1 / 4 + V / 200))
    tau_m = alpha_m * beta_m
    values[3] = (-m + m_inf) / tau_m

    # Expressions for the h gate component
    h_inf = 1.0 / (
        (1 + 15212.593285654404 * math.exp(0.13458950201884254 * V))
        * (1 + 15212.593285654404 * math.exp(0.13458950201884254 * V))
    )
    alpha_h = (
        4.4312679295805147e-07 * math.exp(-0.14705882352941177 * V) if V < -40 else 0
    )
    beta_h = (
        310000 * math.exp(0.3485 * V) + 2.7 * math.exp(0.079 * V)
        if V < -40
        else 0.77 / (0.13 + 0.049758141083938695 * math.exp(-0.0900900900900901 * V))
    )
    tau_h = 1.0 / (alpha_h + beta_h)
    values[4] = (-h + h_inf) / tau_h

    # Expressions for the j gate component
    j_inf = 1.0 / (
        (1 + 15212.593285654404 * math.exp(0.13458950201884254 * V))
        * (1 + 15212.593285654404 * math.exp(0.13458950201884254 * V))
    )
    alpha_j = (
        (37.78 + V)
        * (-25428 * math.exp(0.2444 * V) - 6.948e-06 * math.exp(-0.04391 * V))
        / (1 + 50262745825.95399 * math.exp(0.311 * V))
        if V < -40
        else 0
    )
    beta_j = (
        0.02424
        * math.exp(-0.01052 * V)
        / (1 + 0.003960868339904256 * math.exp(-0.1378 * V))
        if V < -40
        else 0.6 * math.exp(0.057 * V) / (1 + 0.040762203978366204 * math.exp(-0.1 * V))
    )
    tau_j = 1.0 / (alpha_j + beta_j)
    values[5] = (-j + j_inf) / tau_j

    # Expressions for the Sodium background current component
    i_b_Na = g_bna * (-E_Na + V)

    # Expressions for the L_type Ca current component
    V_eff = 0.01 if math.fabs(-15 + V) < 0.01 else -15 + V
    i_CaL = (
        4
        * g_CaL
        * (F * F)
        * (-Ca_o + 0.25 * Ca_ss * math.exp(2 * F * V_eff / (R * T)))
        * V_eff
        * d
        * f
        * f2
        * fCass
        / (R * T * (-1 + math.exp(2 * F * V_eff / (R * T))))
    )

    # Expressions for the d gate component
    d_inf = 1.0 / (1 + 0.34415378686541237 * math.exp(-0.13333333333333333 * V))
    alpha_d = 0.25 + 1.4 / (1 + math.exp(-35 / 13 - V / 13))
    beta_d = 1.4 / (1 + math.exp(1 + V / 5))
    gamma_d = 1.0 / (1 + math.exp(5 / 2 - V / 20))
    tau_d = alpha_d * beta_d + gamma_d
    values[6] = (-d + d_inf) / tau_d

    # Expressions for the f gate component
    f_inf = 1.0 / (1 + math.exp(20 / 7 + V / 7))
    tau_f = (
        20
        + 180 / (1 + math.exp(3 + V / 10))
        + 200 / (1 + math.exp(13 / 10 - V / 10))
        + 1102.5 * math.exp(-((27 + V) * (27 + V)) / 225)
    )
    values[7] = (-f + f_inf) / tau_f

    # Expressions for the F2 gate component
    f2_inf = 0.33 + 0.67 / (1 + math.exp(5 + V / 7))
    tau_f2 = (
        31 / (1 + math.exp(5 / 2 - V / 10))
        + 80 / (1 + math.exp(3 + V / 10))
        + 562 * math.exp(-((27 + V) * (27 + V)) / 240)
    )
    values[8] = (-f2 + f2_inf) / tau_f2

    # Expressions for the FCass gate component
    fCass_inf = 0.4 + 0.6 / (1 + 400.0 * (Ca_ss * Ca_ss))
    tau_fCass = 2 + 80 / (1 + 400.0 * (Ca_ss * Ca_ss))
    values[9] = (-fCass + fCass_inf) / tau_fCass

    # Expressions for the Calcium background current component
    i_b_Ca = g_bca * (-E_Ca + V)

    # Expressions for the Transient outward current component
    i_to = g_to * (-E_K + V) * r * s

    # Expressions for the s gate component
    s_inf = 1.0 / (1 + math.exp(4 + V / 5))
    tau_s = (
        3 + 5 / (1 + math.exp(-4 + V / 5)) + 85 * math.exp(-((45 + V) * (45 + V)) / 320)
    )
    values[10] = (-s + s_inf) / tau_s

    # Expressions for the r gate component
    r_inf = 1.0 / (1 + math.exp(10 / 3 - V / 6))
    tau_r = 0.8 + 9.5 * math.exp(-((40 + V) * (40 + V)) / 1800)
    values[11] = (-r + r_inf) / tau_r

    # Expressions for the Sodium potassium pump current component
    i_NaK = (
        K_o
        * P_NaK
        * Na_i
        / (
            (K_mNa + Na_i)
            * (K_mk + K_o)
            * (
                1
                + 0.0353 * math.exp(-F * V / (R * T))
                + 0.1245 * math.exp(-0.1 * F * V / (R * T))
            )
        )
    )

    # Expressions for the Sodium calcium exchanger current component
    i_NaCa = (
        K_NaCa
        * (
            Ca_o * (Na_i * Na_i * Na_i) * math.exp(F * gamma * V / (R * T))
            - alpha
            * (Na_o * Na_o * Na_o)
            * Ca_i
            * math.exp(F * (-1 + gamma) * V / (R * T))
        )
        / (
            (1 + K_sat * math.exp(F * (-1 + gamma) * V / (R * T)))
            * (Ca_o + Km_Ca)
            * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o))
        )
    )

    # Expressions for the Calcium pump current component
    i_p_Ca = g_pCa * Ca_i / (K_pCa + Ca_i)

    # Expressions for the Potassium pump current component
    i_p_K = (
        g_pK * (-E_K + V) / (1 + 65.40521574193832 * math.exp(-0.16722408026755853 * V))
    )

    # Expressions for the Calcium dynamics component
    i_up = Vmax_up / (1 + (K_up * K_up) / (Ca_i * Ca_i))
    i_leak = V_leak * (-Ca_i + Ca_SR)
    i_xfer = V_xfer * (-Ca_i + Ca_ss)
    kcasr = max_sr - (max_sr - min_sr) / (1 + (EC * EC) / (Ca_SR * Ca_SR))
    Ca_i_bufc = 1.0 / (1 + Buf_c * K_buf_c / ((K_buf_c + Ca_i) * (K_buf_c + Ca_i)))
    Ca_sr_bufsr = 1.0 / (
        1 + Buf_sr * K_buf_sr / ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR))
    )
    Ca_ss_bufss = 1.0 / (
        1 + Buf_ss * K_buf_ss / ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss))
    )
    values[12] = (
        V_sr * (-i_up + i_leak) / V_c
        - Cm * (-2 * i_NaCa + i_b_Ca + i_p_Ca) / (2 * F * V_c)
        + i_xfer
    ) * Ca_i_bufc
    k1 = k1_prime / kcasr
    k2 = k2_prime * kcasr
    O = (Ca_ss * Ca_ss) * R_prime * k1 / (k3 + (Ca_ss * Ca_ss) * k1)  # noqa: E741
    values[13] = k4 * (1 - R_prime) - Ca_ss * R_prime * k2
    i_rel = V_rel * (-Ca_ss + Ca_SR) * O
    values[14] = (-i_leak - i_rel + i_up) * Ca_sr_bufsr
    values[15] = (
        V_sr * i_rel / V_ss - V_c * i_xfer / V_ss - Cm * i_CaL / (2 * F * V_ss)
    ) * Ca_ss_bufss

    # Expressions for the Sodium dynamics component
    values[16] = Cm * (-i_Na - i_b_Na - 3 * i_NaCa - 3 * i_NaK) / (F * V_c)

    # Expressions for the Membrane component
    i_Stim = (
        -stim_amplitude
        if t - stim_period * math.floor(t / stim_period) <= stim_duration + stim_start
        and t - stim_period * math.floor(t / stim_period) >= stim_start
        else 0
    )
    values[17] = (
        -i_CaL
        - i_K1
        - i_Kr
        - i_Ks
        - i_Na
        - i_NaCa
        - i_NaK
        - i_Stim
        - i_b_Ca
        - i_b_Na
        - i_p_Ca
        - i_p_K
        - i_to
    )

    # Expressions for the Potassium dynamics component
    values[18] = (
        Cm * (-i_K1 - i_Kr - i_Ks - i_Stim - i_p_K - i_to + 2 * i_NaK) / (F * V_c)
    )

    # Return results
    return values


selected_methods = [
    "goss (ExplicitEuler)",
    "goss (RK2)",
    "goss (RK4)",
    "goss (RL1)",
    "goss (RL2)",
    "goss (GRL1)",
    "goss (GRL2)",
    "goss (ImplicitEuler)",
    "goss (ThetaSolver)",
    "goss (RKF32)",
    "goss (ESDIRK23a)",
    # "scipy (RK45)",
    # "scipy+numba (RK45)",
    # "scipy (RK23)",
    # "scipy+numba (RK23)",
    # "scipy (DOP853)",
    # "scipy+numba (DOP853)",
    # "scipy (Radau)",
    # "scipy+numba (Radau)",
    # "scipy (BDF)",
    # "scipy+numba (BDF)",
    # "scipy (LSODA)",
    # "scipy+numba (LSODA)",
]

internal_time_steps = benchmark.default_internal_time_steps()
internal_time_steps["goss (ExplicitEuler)"] = 0.001
internal_time_steps["goss (RK2)"] = 0.01
internal_time_steps["goss (RK4)"] = 0.01
internal_time_steps["goss (RL1)"] = 0.01
internal_time_steps["goss (RL2)"] = 0.01
internal_time_steps["goss (GRL1)"] = 0.01
internal_time_steps["goss (GRL2)"] = 0.1
internal_time_steps["goss (ImplicitEuler)"] = 0.01
internal_time_steps["goss (ThetaSolver)"] = 0.1
# internal_time_steps["goss (RKF32)"] = 0.001
# internal_time_steps["goss (ESDIRK23a)"] = 0.001


def main():

    init_parameters = np.array(
        [
            0.03,
            5.405,
            0.153,
            0.098,
            14.838,
            0.00029,
            3.98e-05,
            0.000592,
            0.294,
            40,
            1,
            2.724,
            1000,
            0.1,
            1.38,
            87.5,
            2.5,
            0.35,
            0.0005,
            0.1238,
            0.0146,
            0.2,
            10,
            0.4,
            2,
            1.5,
            0.001,
            0.3,
            0.00025,
            0.00025,
            0.00036,
            0.102,
            0.001094,
            5.468e-05,
            0.0038,
            0.006375,
            0.15,
            0.045,
            0.06,
            0.005,
            2.5,
            1.0,
            140,
            0.185,
            96485.3415,
            8314.472,
            310,
            0.016404,
            52,
            1,
            1000,
            10,
            5.4,
        ],
        dtype=np.float_,
    )
    benchmark.main(
        rhs,
        "tentusscher_panfilov_2006_M_cell.ode",
        t=np.linspace(0, 1000, num=1001),
        args=(init_parameters,),
        run_plot=True,
        run_timings=False,
        recompute=True,
        selected_methods=selected_methods,
        selected_states=None,
        # selected_states=[
        #     12,
        #     14,
        #     17,
        #     18,
        # ],
        internal_time_step=internal_time_steps,
        number=2,
        repeat=4,
    )


if __name__ == "__main__":
    main()
