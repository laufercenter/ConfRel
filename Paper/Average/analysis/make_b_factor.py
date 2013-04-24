#!/usr/bin/env python
# encoding: utf-8

import prody
import numpy
import math


def smooth_b_factors(protein, cutoff):
    return protein

    energies = []

    for i, atom_i in enumerate(protein):
        w_total = 0.
        energy_total = 0.
        for j, atom_j in enumerate(protein):
            r = numpy.linalg.norm(atom_i.getCoords() - atom_j.getCoords())
            w = math.exp(-r ** 2 / (0.5 * cutoff ** 2))
            energy_total += w * atom_j.getBeta()
            w_total += w
        energies.append(energy_total / w_total)

    protein.setBetas(energies)
    return protein


def set_betas(energies, protein):
    for i, atom in enumerate(protein):
        res_index = atom.getResnum() - 1
        atom.setBeta(energies[res_index])
    return protein


def get_per_residue(filename):
    energies = numpy.loadtxt(filename)[:, 1:]
    energies[:, -1] = 0.0072 * energies[:, -1]
    return numpy.sum(energies, axis=1)

energies_ga_a = -1 * numpy.loadtxt('confinement_Alpha_GA95')[:, 1] + 1 * numpy.loadtxt('enthalpy_spring_Alpha_GA95')[:, 1]
energies_ga_a += get_per_residue('residue_Alpha_GA95')
alpha_ga = prody.parsePDB('../Alpha/GA95/A95_min2.pdb', subset='CA')
alpha_ga = set_betas(energies_ga_a, alpha_ga)
alpha_ga = smooth_b_factors(alpha_ga, 5.0)

energies_ga_b = -1 * numpy.loadtxt('confinement_Beta_GA95')[:, 1] + 1 * numpy.loadtxt('enthalpy_spring_Beta_GA95')[:, 1]
energies_ga_b += 1 * get_per_residue('residue_Beta_GA95')
beta_ga = prody.parsePDB('../Beta/GA95/A95_min.pdb', subset='CA')
beta_ga = set_betas(energies_ga_b, beta_ga)
beta_ga = smooth_b_factors(beta_ga, 5.0)

energies_gb_a = -1 * numpy.loadtxt('confinement_Alpha_GB95')[:, 1] + 1 * numpy.loadtxt('enthalpy_spring_Alpha_GB95')[:, 1]
energies_gb_a += 1 * get_per_residue('residue_Alpha_GB95')
alpha_gb = prody.parsePDB('../Alpha/GB95/Alpha_GB95.pdb', subset='CA')
alpha_gb = set_betas(energies_gb_a, alpha_gb)
alpha_gb = smooth_b_factors(alpha_gb, 5.0)

energies_gb_b = -1 * numpy.loadtxt('confinement_Beta_GB95')[:, 1] + 1 * numpy.loadtxt('enthalpy_spring_Beta_GB95')[:, 1]
energies_gb_b += 1 * get_per_residue('residue_Beta_GB95')
beta_gb = prody.parsePDB('../Beta/GB95/B95_min.pdb', subset='CA')
beta_gb = set_betas(energies_gb_b, beta_gb)
beta_gb = smooth_b_factors(beta_gb, 5.0)

betas_ga = beta_ga.getBetas() - alpha_ga.getBetas()
betas_gb = beta_gb.getBetas() - alpha_gb.getBetas()
#betas_ga = betas_ga - betas_gb
#betas_gb = betas_ga.copy()

#cutoff = 2.5
#betas_ga[ numpy.abs(betas_ga) < cutoff ] = 0.
#betas_ga[ betas_ga < -cutoff ] = -1.
#betas_ga[ betas_ga > cutoff ] = 1.
#betas_gb[ numpy.abs(betas_gb) < cutoff ] = 0.
#betas_gb[ betas_gb < -cutoff ] = -1.
#betas_gb[ betas_gb > cutoff ] = 1.


alpha_ga = prody.parsePDB('../Alpha/GA95/A95_min2.pdb')
beta_ga = prody.parsePDB('../Beta/GA95/A95_min.pdb')
alpha_gb = prody.parsePDB('../Alpha/GB95/Alpha_GB95.pdb')
beta_gb = prody.parsePDB('../Beta/GB95/B95_min.pdb')

alpha_ga = set_betas(betas_ga, alpha_ga)
beta_ga = set_betas(betas_ga, beta_ga)
alpha_gb = set_betas(betas_gb, alpha_gb)
beta_gb = set_betas(betas_gb, beta_gb)

prody.writePDB('alpha_ga.pdb', alpha_ga)
prody.writePDB('beta_ga.pdb', beta_ga)
prody.writePDB('alpha_gb.pdb', alpha_gb)
prody.writePDB('beta_gb.pdb', beta_gb)
