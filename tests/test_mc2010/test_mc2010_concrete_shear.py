"""Test for the function of _concrete_shear."""

import math

import pytest

from structuralcodes.codes.mc2010 import _concrete_shear


def create_load_dict(
    Med: float, Ved: float, Ned: float, delta_e: float
) -> dict:
    """Returns dictionary assosiated with loads."""
    return {'Med': Med, 'Ved': Ved, 'Ned': Ned, 'delta_e': delta_e}


@pytest.mark.parametrize(
    'E_s, As, z, loads, expected',
    [
        (
            200000,
            1000,
            160,
            create_load_dict(50000000, 10000, 2000, 50),
            8.1e-4,
        ),
        (
            210000,
            1000,
            160,
            create_load_dict(50000000, 10000, 2000, 50),
            7.7e-4,
        ),
        (
            210000,
            5000,
            160,
            create_load_dict(50000000, 10000, 2000, 50),
            1.5e-4,
        ),
        (
            210000,
            2000,
            160,
            create_load_dict(50000000, 10000, 2000, 50),
            3.9e-4,
        ),
        (
            210000,
            2000,
            160,
            create_load_dict(40000000, 20000, 2000, 50),
            3.2e-4,
        ),
        (
            210000,
            2000,
            160,
            create_load_dict(40000000, 20000, 1000, 50),
            3.2e-4,
        ),
        (
            210000,
            2000,
            140,
            create_load_dict(40000000, 20000, 1000, 50),
            3.64965e-4,
        ),
        (
            210000,
            2000,
            180,
            create_load_dict(40000000, 20000, 1000, 50),
            2.9e-4,
        ),
    ],
)
def test_epsilon_x(E_s, As, z, loads, expected):
    """Test the epsilon_x function."""
    assert math.isclose(
        _concrete_shear.epsilon_x(E_s, As, z, loads), expected, rel_tol=0.05
    )


@pytest.mark.parametrize(
    'approx_lvl, fck, bw, theta, z, E_s, As, loads, alpha, gamma_c, expected',
    [
        (
            1,
            30,
            50,
            20,
            200,
            210000,
            1000,
            create_load_dict(200e6, 50e3, 10e3, 50),
            20,
            1.5,
            70707,
        ),
        (
            2,
            30,
            50,
            20,
            200,
            210000,
            1000,
            create_load_dict(200e6, 50e3, 10e3, 50),
            20,
            1.5,
            39997,
        ),
        (
            2,
            30,
            50,
            20,
            200,
            210000,
            1000,
            create_load_dict(50e6, 10e3, 10e3, 50),
            20,
            1.5,
            55179.55,
        ),
        (
            2,
            30,
            50,
            45,
            200,
            210000,
            1000,
            create_load_dict(0, 0, 0, 50),
            20,
            1.5,
            243586,
        ),
        (
            2,
            30,
            50,
            45,
            200,
            210000,
            1000,
            create_load_dict(0, 0, 0, 50),
            45,
            1.5,
            130000,
        ),
        (
            3,
            30,
            50,
            20,
            200,
            210000,
            1000,
            create_load_dict(50e6, 10e3, 10e3, 50),
            20,
            1.5,
            102995,
        ),
    ],
)
def test_vrd_max(
    approx_lvl, fck, bw, theta, z, E_s, As, loads, alpha, gamma_c, expected
):
    """Test the v_rd_max function."""
    assert math.isclose(
        _concrete_shear.v_rd_max(
            approx_lvl, fck, bw, theta, z, E_s, As, loads, alpha, gamma_c
        ),
        expected,
        rel_tol=0.5,
    )


@pytest.mark.parametrize(
    'approx_lvl, fck, z, bw, dg, E_s, As, loads, alpha, gamma_c, expected',
    [
        (
            1,
            35,
            180,
            300,
            0,
            0,
            0,
            create_load_dict(0, 0, 0, 0),
            0,
            1.5,
            31294,
        ),
        (
            1,
            35,
            200,
            300,
            0,
            0,
            0,
            create_load_dict(0, 0, 0, 0),
            0,
            1.5,
            34077,
        ),
        (
            2,
            35,
            140,
            300,
            16,
            21e4,
            2000,
            create_load_dict(40e6, 2e4, 1000, 50),
            0,
            1.5,
            48828,
        ),
        (
            2,
            35,
            140,
            300,
            32,
            21e4,
            2000,
            create_load_dict(40e6, 2e4, 1000, 50),
            0,
            1.5,
            50375,
        ),
        (
            3,
            35,
            200,
            300,
            32,
            21e4,
            2000,
            create_load_dict(40e6, 2e4, 1000, 50),
            1.5,
            1.5,
            67566,
        ),
        (
            3,
            35,
            200,
            300,
            32,
            21e4,
            2000,
            create_load_dict(40e6, 20e6, 1000, 50),
            1.5,
            1.5,
            0,
        ),
    ],
)
def test_v_rdc(
    approx_lvl, fck, z, bw, dg, E_s, As, loads, alpha, gamma_c, expected
):
    """Test the v_rdc function."""
    assert math.isclose(
        _concrete_shear.v_rdc(
            approx_lvl, fck, z, bw, dg, E_s, As, loads, alpha, gamma_c
        ),
        expected,
        rel_tol=0.001,
    )


@pytest.mark.parametrize(
    'asw, sw, z, fywd, theta, alpha, expected',
    [
        (1600, 50, 200, 500, 25, 30, 5383763),
        (2000, 50, 200, 500, 25, 30, 6729704),
        (1600, 50, 200, 500, 25, 30, 5383763),
        (1600, 100, 200, 500, 25, 30, 2691881),
        (1600, 50, 200, 500, 25, 30, 5383763),
        (1600, 50, 200, 500, 22, 30, 5842872),
        (1600, 50, 200, 500, 25, 25, 5034721),
    ],
)
def test_v_rds(asw, sw, z, fywd, theta, alpha, expected):
    """Test the v_rds function."""
    assert math.isclose(
        _concrete_shear.v_rds(asw, sw, z, fywd, theta, alpha),
        expected,
        rel_tol=0.005,
    )


@pytest.mark.parametrize(
    (
        'approx_lvl_h, f_ctd, i_c, s_c, b_w, sigma_cp, l_x, l_bd0, S_cy, b_wy,'
        ' y, y_c, A_c, A_cy, y_pt, f_p_lx, f_p_lx_dx, expected'
    ),
    [
        (
            1,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            839136,
        ),
        (
            1,
            3.5,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            976183,
        ),
        (
            1,
            2.6,
            5e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            699280,
        ),
        (
            1,
            2.6,
            6e8,
            5e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            1006963,
        ),
        (
            1,
            2.6,
            6e8,
            6e5,
            40,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            671309,
        ),
        (
            1,
            2.6,
            6e8,
            6e5,
            50,
            180,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            918050,
        ),
        (
            1,
            2.6,
            6e8,
            6e5,
            50,
            150,
            35,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            785800,
        ),
        (
            1,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            25,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            918050,
        ),
        (
            1,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            2e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            839136,
        ),
        (
            1,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            2e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            839136,
        ),
        (
            1,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            80,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            839136,
        ),
        (
            1,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            180,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            839136,
        ),
        (
            1,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            1800,
            1000,
            80,
            1000e3,
            200e3,
            839136,
        ),
        (
            1,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1200,
            80,
            1000e3,
            200e3,
            839136,
        ),
        (
            1,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            60,
            1000e3,
            200e3,
            839136,
        ),
        (
            1,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            800e3,
            200e3,
            839136,
        ),
        (
            1,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            250e3,
            839136,
        ),
        (
            2,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            160002777,
        ),
        (
            2,
            3.5,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            214002777,
        ),
        (
            2,
            2.6,
            5e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            137336111,
        ),
        (
            2,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            160002777,
        ),
        (
            2,
            2.6,
            6e8,
            6e5,
            40,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            160002777,
        ),
        (
            2,
            2.6,
            6e8,
            6e5,
            50,
            180,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            160002777,
        ),
        (
            2,
            2.6,
            6e8,
            6e5,
            50,
            150,
            35,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            160002430,
        ),
        (
            2,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            25,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            160003333,
        ),
        (
            2,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            2e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            228004166,
        ),
        (
            2,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            2e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            108001851,
        ),
        (
            2,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            80,
            200,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            160003333,
        ),
        (
            2,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            180,
            2000,
            1000,
            80,
            1000e3,
            200e3,
            156002222,
        ),
        (
            2,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            1800,
            1000,
            80,
            1000e3,
            200e3,
            157780864,
        ),
        (
            2,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1200,
            80,
            1000e3,
            200e3,
            156002777,
        ),
        (
            2,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            60,
            1000e3,
            200e3,
            164002777,
        ),
        (
            2,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            800e3,
            200e3,
            160002222,
        ),
        (
            2,
            2.6,
            6e8,
            6e5,
            50,
            150,
            40,
            30,
            3e6,
            3e5,
            100,
            200,
            2000,
            1000,
            80,
            1000e3,
            250e3,
            161002777,
        ),
    ],
)
def test_v_rd_ct(
    approx_lvl_h,
    f_ctd,
    i_c,
    s_c,
    b_w,
    sigma_cp,
    l_x,
    l_bd0,
    S_cy,
    b_wy,
    y,
    y_c,
    A_c,
    A_cy,
    y_pt,
    f_p_lx,
    f_p_lx_dx,
    expected,
):
    """Test the v_rd_ct function."""
    assert math.isclose(
        _concrete_shear.v_rd_ct(
            approx_lvl_h,
            f_ctd,
            i_c,
            s_c,
            b_w,
            sigma_cp,
            l_x,
            l_bd0,
            S_cy,
            b_wy,
            y,
            y_c,
            A_c,
            A_cy,
            y_pt,
            f_p_lx,
            f_p_lx_dx,
        ),
        expected,
        rel_tol=0.001,
    )


@pytest.mark.parametrize(
    (
        'approx_lvl, with_shear_reinforcment, fck, z, bw, dg, E_s, As, loads,'
        ' asw, sw, f_ywk, theta, alpha, gamma_c, gamma_s, expected'
    ),
    [
        (
            1,
            False,
            35,
            180,
            200,
            16,
            200000,
            2000,
            create_load_dict(0, 2000, 0, 20),
            0,
            0,
            500,
            40,
            90,
            1.5,
            1.15,
            20863,
        ),
        (
            2,
            False,
            35,
            180,
            200,
            16,
            200000,
            2000,
            create_load_dict(0, 2000, 0, 20),
            0,
            0,
            500,
            40,
            90,
            1.5,
            1.15,
            62336,
        ),
        (
            1,
            True,
            35,
            180,
            200,
            16,
            200000,
            2000,
            create_load_dict(0, 2000, 0, 20),
            500,
            200,
            500,
            40,
            90,
            1.5,
            1.15,
            216096,
        ),
        (
            2,
            True,
            35,
            180,
            200,
            16,
            200000,
            2000,
            create_load_dict(0, 2000, 0, 20),
            500,
            200,
            500,
            40,
            90,
            1.5,
            1.15,
            233170,
        ),
        (
            3,
            True,
            35,
            180,
            200,
            16,
            200000,
            2000,
            create_load_dict(0, 2000, 0, 20),
            500,
            200,
            500,
            40,
            90,
            1.5,
            1.15,
            233170,
        ),
    ],
)
def test_v_rd(
    approx_lvl,
    with_shear_reinforcment,
    fck,
    z,
    bw,
    dg,
    E_s,
    As,
    loads,
    asw,
    sw,
    f_ywk,
    theta,
    alpha,
    gamma_c,
    gamma_s,
    expected,
):
    """Test the tau_edi function."""
    assert math.isclose(
        _concrete_shear.v_rd(
            approx_lvl,
            with_shear_reinforcment,
            fck,
            z,
            bw,
            dg,
            E_s,
            As,
            loads,
            asw,
            sw,
            f_ywk,
            theta,
            alpha,
            gamma_c,
            gamma_s,
        ),
        expected,
        rel_tol=0.001,
    )
