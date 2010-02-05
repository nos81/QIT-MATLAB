% Guided tour to the quantum information toolkit.

% Ville Bergholm 2009-2010

disp('This is the guided tour for the Quantum Information Toolkit.')
disp('Between examples, press any key to proceed to the next one.');
pause

teleportation;
pause

superdense_coding;
pause

adiabatic_qc(5, 3);
pause

phase_estimation(5, rand_U(4));
title('Phase estimation, eigenstate')
phase_estimation(5, rand_U(4), state(0, [4]));
title('Phase estimation, random state')
pause

nmr_sequences;
pause

quantum_channels(0.3)
pause

grover_search(6);
pause

bb84(40);
pause

markov_decoherence(7e-10, 1e-9);
pause

qubit_and_resonator();
