# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2018 Allan Leal
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.

from reaktoro import *

# Auxiliary constants
mNaCl = 5.0  # the molality of NaCl used below

# Define the chemical system
db = Database('supcrt98.xml')

editor = ChemicalEditor(db)

editor.addAqueousPhase(['H2O(l)', 'H+', 'OH-', 'Cl-', 'Na+']) \
    .setChemicalModelPitzerHMW()

system = ChemicalSystem(editor)

# Set a equilibrium problem for pure water
problem = EquilibriumProblem(system)
problem.setTemperature(25, "celsius")
problem.setPressure(1, "bar")
problem.add("H2O", 1, "kg")

# Compute the equilibrium state of pure water
state_pw = equilibrate(problem)

# Get the chemical potential of water solvent in pure water solution
properties = state_pw.properties()

u0 = properties.standardPartialMolarGibbsEnergies().val
u = properties.chemicalPotentials().val
iH2O = system.indexSpecies("H2O(l)")

u0H2O, uH2O = u0[iH2O], u[iH2O]

# Perform an equilibrium calculation in which solvent water
# in a saline solution has a fixed chemical potential of that
# of pure water, and pressure is unknown.
problem = EquilibriumInverseProblem(system)
problem.add("H2O", 1, "kg")
problem.add("NaCl", 5, "mol")
problem.setTemperature(25, "celsius")
problem.setPressure(1, "bar")
problem.unknownPressure()
problem.fixSpeciesChemicalPotential("H2O(l)", uH2O)

# Print the iterations of the inverse equilibrium calculation
options = EquilibriumOptions()
options.nonlinear.output.active = True

# Solve the inverse chemical equilibrium calculation
state = equilibrate(problem, options)

# Output the chemical state to a file
state.output('demo-equilibrium-fixed-chemical-potential-water-unknown-pressure.txt')

# Compute the osmotic coefficient of solvent water
properties = state.properties()
ln_aw = properties.lnActivities().val[iH2O]
phi = -ln_aw / (2*mNaCl*waterMolarMass)

print('Osmotic coefficient of water: {}'.format(phi))
