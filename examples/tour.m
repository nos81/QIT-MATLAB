% Guided tour to the quantum information toolkit.

% Ville Bergholm 2009-2011

disp('This is the guided tour for the Quantum Information Toolkit.')
disp('Between examples, press any key to proceed to the next one.');
pause

teleportation(2);
pause

superdense_coding(2);
pause

adiabatic_qc_3sat(5, 25);
pause

phase_estimation_precision(5, rand_U(4));
title('Phase estimation, eigenstate')
phase_estimation_precision(5, rand_U(4), state(0, [4]));
title('Phase estimation, random state')
pause

nmr_sequences;
pause

quantum_channels(0.3)
pause

bernstein_vazirani(6, true);
bernstein_vazirani(6, false);
pause

grover_search(6);
pause

shor_factorization(9);
shor_factorization(311*269, true);
pause

bb84(40);
pause

markov_decoherence(7e-10, 1e-9);
pause

qubit_and_resonator();
pause

qft_circuit([2, 3, 2]);
pause

werner_states(2);
werner_states(3);
