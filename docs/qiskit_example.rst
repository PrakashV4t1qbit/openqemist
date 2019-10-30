Using Qiskit with QEMIST
=======================

VQE Example
-----------

This tutorial assumes that the user has correctly set up and configured
the OpenQEMIST package. The Variational Quantum Eigensolver
(VQE) :math:`^{1,2}` is a hybrid quantum-classical algorithm for
simulating quantum systems. We here focus on VQE within the context of
solving the molecular electronic structure problem for the ground-state
energy of a molecular system. In VQE, we first prepare the trial
wavefunction (quantum state)
:math:`\vert \Psi(\vec{\theta}) \rangle = U(\vec{\theta}) \vert 0 \rangle`
based on an ansatz that depends on :math:`m` parameters defining
:math:`\vec{\theta}=(\theta_1, \theta_2, \ldots, \theta_m)`. The
expectation value of the Hamiltonian (:math:`\hat{H}`),
:math:`\langle \Psi(\vec{\theta}) \vert \hat{H} \vert \Psi(\vec{\theta}) \rangle`,
will then be simulated.

The expectation value can be minimized based on the variational
principle,

.. raw:: latex

   \begin{equation}
   E = \min_{\vec{\theta}} \frac{\langle \Psi(\vec{\theta}) \vert \hat{H} \vert \Psi(\vec{\theta}) \rangle}{\langle \Psi(\vec{\theta}) \vert \Psi(\vec{\theta}) \rangle} \geq E_{\text{gs}}\nonumber
   \end{equation}

which ensures that the energy computed will be an upper bound to the
true ground-state energy :math:`E_{\text{gs}}`. This allows us using
classical minimizers to find optimal parameters :math:`\vec{\theta}` for
the ground-state energy :math:`E_{\text{gs}}`.

VQE can be performed using OpenQEMIST in conjuction with Qiskit for
calculating the ground state energy of a molecular system. Unitary
coupled-cluster ansatz can be used to prepare the trial wavefunction
:math:`\vert \Psi(\vec{\theta}) \rangle`. In this notebook, we will show
you an example using a small molecule, the hydrogen molecule
(H :math:`_\text{2}`), for a simulation using VQE.

.. code:: ipython3

    from openqemist.electronic_structure_solvers import VQESolver
    from openqemist.quantum_solvers import QiskitParametricSolver
    
    from pyscf import gto
    
    # Build the molecule
    H2 = [['H', [0.0,   0.0,   0.0]], ['H', [0.0,   0.0,   0.74137727]]]
    mol = gto.Mole()
    mol.atom = H2
    mol.basis = "minao"
    mol.charge = 0
    mol.spin = 0
    mol.build()
    
    # Configure the solver object
    solver = VQESolver()
    solver.hardware_backend_type = QiskitParametricSolver
    solver.ansatz_type = QiskitParametricSolver.Ansatze.UCCSD

We can now simulate the molecule and get its energy.

.. code:: ipython3

    energy = solver.simulate(mol)


.. parsed-literal::

    VQE : initial variational parameters: 
     [0. 0. 0.] 
    
.. parsed-literal::

    Optimization terminated successfully.    (Exit mode 0)
                Current function value: -1.1062746291459193
                Iterations: 2
                Function evaluations: 11
                Gradient evaluations: 2
    
    		Optimal UCCSD Singlet Energy: -1.1062746291459193
    		Optimal UCCSD Singlet Amplitudes: [-2.38337414e-06 -2.38337414e-06 -1.11078091e-01]
    		Number of Function Evaluations :  11


It is possible to use different initial parameters for the optimization:

.. code:: ipython3

    # Using custom initial parameters
    # Getting the dimension of the initial parameters vector
    num_var_params = solver.hardware_backend.amplitude_dimension
    # Set the intial parameters for the solver
    solver.initial_var_params = [0.1 for i in range(num_var_params)]
    
    solver.simulate(mol)


.. parsed-literal::

    VQE : initial variational parameters: 
     [0.1, 0.1, 0.1] 

.. parsed-literal::

    Optimization terminated successfully.    (Exit mode 0)
                Current function value: -1.1062743622818396
                Iterations: 4
                Function evaluations: 21
                Gradient evaluations: 4
    
    		Optimal UCCSD Singlet Energy: -1.1062743622818396
    		Optimal UCCSD Singlet Amplitudes: [ 0.00031692  0.00031692 -0.11127638]
    		Number of Function Evaluations :  21

.. parsed-literal::

    -1.1062743622818396



For advanced usage, including configuring the optimization algorithm,
please refer to the reference manual or the other examples included with
OpenQEMIST.

DMET Example
------------

At the current early stage of quantum hardware, the available
computational resource is yet very limited. Thus, it is still
challenging to perform accurate electronic structure calculations on
actual quantum hardware. Simulation on classical computer requires large
computational cost as well. Therefore, we need to reduce the problem
size while maintaining the accuracy of electronic structure calculation
to solve a problem for small sized molecules to perform quantum
simulations.

Density Matrix Embedding Theory (DMET) :math:`^{3,4}` is a powerful
problem decomposition technique to reduce the problem size, while
maintaining the accuracy of the electronic structure calculation. The
DMET method decomposes a molecule into fragments, and each fragment is
treated as an open quantum system that is entangled with each of the
other fragments, all taken together to be that fragment’s surrounding
environment (or “bath”). VQE algorithm can be used with DMET using
OpenQEMIST in conjuction with Qiskit.

In this notebook, we will show you an example of H\ :math:`_\text{4}`
molecule for DMET simulation using VQE as an electronic structure
solver.

.. code:: ipython3

    from openqemist.problem_decomposition import DMETProblemDecomposition
    from openqemist.problem_decomposition.electron_localization import meta_lowdin_localization
    
    
    H4 = [['H', [0.7071067811865476,   0.0,                 0.0]],
          ['H', [0.0,                  0.7071067811865476,  0.0]],
          ['H', [-1.0071067811865476,  0.0,                 0.0]],
          ['H', [0.0,                 -1.0071067811865476,  0.0]]]
    
    mol = gto.Mole()
    mol.atom = H4
    mol.basis = "minao"
    mol.charge = 0
    mol.spin = 0
    mol.build()
    
    dmet = DMETProblemDecomposition()
    dmet.electron_localization_method = meta_lowdin_localization
    # Set the DMET object to use the solver that we configured above
    dmet.electronic_structure_solver = solver
    energy_vqe = dmet.simulate(mol, [2,2])
    
    print("DMET energy is: ", energy)

.. parsed-literal::

    DMET energy is:  -1.1062746291459193


.. parsed-literal::

      warnings.warn(msg, RuntimeWarning)


References
----------

1. Alberto Peruzzo, Jarrod McClean, Peter Shadbolt, Man-Hong Yung,
   Xiao-Qi Zhou, Peter J. Love, Alán Aspuru-Guzik, and Jeremy L.
   O’Brien, “A variational eigenvalue solver on a photonic quantum
   processor”, Nat. Commun., 5, 4213 (2014).
2. Jarrod R. McClean, Jonathan Romero, Ryan Babbush, and Alán
   Aspuru-Guzik, “The theory of variational hybrid quantum-classical
   algorithms”, New J. Phys., 18, 023023 (2016).
3. Gerald Knizia and Garnet K.-L. Chan, “Density Matrix Embedding: A
   Simple Alternative to Dynamical Mean-Field Theory”, Phys. Rev. Lett.,
   109, 186404 (2012).
4. Sebastian Wouters, Carlos A. Jiménez-Hoyos, Qiming Sun, and Garnet
   K.-L. Chan, “A Practical Guide to Density Matrix Embedding Theory in
   Quantum Chemistry”, J. Chem. Theory Comput., 12, pp. 2706–2719
   (2016).
