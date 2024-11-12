1) Open Find_disturbance_ellipsoid_for_Indirect.m and run. This returns the matrix Y in (16) in the paper (line 65);
2) Open Return_K_and_Phi_by_solve_prob_22_indirect_method.m and enter the obtained Y in line 10 and run. This returns K_ind, \Phi and runtime of the simulation.
3) Run Lemma3.m In line 27 change the feedback gain Kind with the output of the previous step
4) Store Cu (line 45)
5) Open Solve_prob_in_eq7_with_indirect.m and in lines 23, 24, 25 enter the Kind, \Phi, and Cu, respectively from the previous steps and run. (requires YALMIP and MOSEK solver, and MPT3 for plots)
6) The last returns the satisfaction frequencies verified over 10000 disturbance sequence samples (change as needed in line 59) and the solve time of MOSEK solver. It also prints Fig. 1 (right).
