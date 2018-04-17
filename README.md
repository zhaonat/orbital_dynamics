# Nonlinear Orbital Dynamics
Simulating complex Newtonian dynamics with the Runge-Kutta algorithm. The Runge Kutta algorithm reproduces the trajectory of a particle in a given potential with 4th order accuracy. In order to write down the algorithm, one needs the explicit Hamiltonian of the system (which must be differentiable so you can derive the equations of motions).

# Runge-Kutta
folder: runge_kutta_integrator
these are a series of supporting matlab functions which implements the RK algorithm for the Henon-Heiles Hamiltonian.
the Hamiltonian has to be modified since the RK algorithm differentiates the Hamiltonian in order to perform the updates
(need symbolic differentiator like Mathematica)

# Henon-Heiles
folder: orbital_dynamics_Henon_Heiles
main programs in orbital dynamics investigates the runge-kutta algorithm to simulate nonlinear Newtonian dynamics of a particle with the Henon-Heiles Hamiltonian. You should be able to execute any of these files when you download the repo

# Examples
you can find some simpler systems such as the SHO or a 1D projectile simulated with the Runge-Kutta algorithm

# Example Phase Diagram
The Henon-Heiles Hamiltonian is one of the first Hamiltonians studied which exhibits chaos. In short, the notion of chaos in classical dynamics describes the fact that a small perturbation of the initial conditions can lead to a large change in the trajectory of the particle after a sufficiently long period of time.
<br>
![](example_figs/E_0.083_phase_diagram.png?raw=true)

# Example Trajectory
