% Guided tour to the quantum information toolkit.

% Ville Bergholm 2009

disp('This is the guided tour for the Quantum Information Toolkit.')
disp('Between examples, press any key to proceed to the next one.');
pause

teleportation;
pause

adiabatic_qc(5, 3);
pause

phase_estimation(4, rand_U(4));
pause

nmr_sequences;
pause

quantum_channels(0.3)
pause

grover_search(6);
pause

bb84(40);
pause