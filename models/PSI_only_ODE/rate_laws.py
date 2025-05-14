import numpy as np


def moiety_1(concentration, total):
    return total - concentration

def moiety_2(c1, c2, total): # first the compounds are passed and then the parameters
    c3 = total - c1 - c2
    return c3

def moiety_3(c1, c2, c3, total): # first the compounds are passed and then the parameters
    c4 = total - c1 - c2 - c3
    return c4

def mass_action_1s(s1, kf):
    return kf * s1

def mass_action_2s(s1, s2, kf):
    return kf * s1 * s2

def mass_action_22_rev(s1, s2, p1, p2, kf, keq): # reverse reaction
    forward = kf * s1 * s2
    reverse = kf/keq * p1 * p2
    return forward - reverse


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


def normalize_concentration(concentration, total):
    return concentration / total

def normalize_2_concentrations(c1, c2, total):
    return (c1+c2) / total