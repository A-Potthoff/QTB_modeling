import numpy as np


def moiety_1(concentration, total):
    return total - concentration


def mass_action_1s(s1, kf):
    return kf * s1


def mass_action_2s(s1, s2, kf):
    return kf * s1 * s2


def proportional(base, factor):
    return base * factor


def michaelis_menten(s, vmax, km):
    return s * vmax / (s + km)


def rapid_eq_1_1(s1, p1, k, q):
    return k * (s1 - p1 / q)


def rapid_eq_2_1(s1, s2, p1, k, q):
    return k * (s1 * s2 - (p1 / q))


def rapid_eq_2_2(s1, s2, p1, p2, k, q):
    return k * (s1 * s2 - (p1 * p2 / q))


def rapid_eq_3_3(s1, s2, s3, p1, p2, p3, k, q):
    return k * (s1 * s2 * s3 - (p1 * p2 * p3 / q))


def vPS1(A, ps2cs, pfd):
    """reaction rate constant for open PSI"""
    return (1 - ps2cs) * pfd * A


def normalize_concentration(concentration, total):
    return concentration / total
        
