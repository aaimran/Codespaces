# README #

WaveQLab3D is a code for 3D seismic wave propagation and earthquake rupture dynamics. It solves the elastic wave equation in curvilinear coordinates (i.e., complex geometries) with a possibly nonplanar frictional fault interface. The current version supports off-fault viscoplasticity, spatially variable elastic properties, and several friction laws (including rate-and-state and slip-weakening). The code is under development and is available under the MIT license. Authors include Kenneth Duru, Sam Bydlon, Eric Dunham, and Kyle Withers with parallelization by Hari Radhakrishnan.

Supported attenuation response options currently include `anelastic`, `anelastic-Q`, `anelastic-Q8`, `anelastic-Qf`, `constant-Q-4M`, `constant-Q-8M`, `frequency-Q-4M`, and `frequency-Q-8M`.
