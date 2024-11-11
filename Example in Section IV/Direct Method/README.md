1) Open the jupyter notebook GA_direct_example_section_IV and run it. It will return the feedback gain. 
2) Run Lemma3.m in matlab. In line 28 change the feedback gain Kdir with the output of the previous step(GA)
3) Store [Ce Cu] (line 46)
4) Open Solve_prob_in_eq7_with_direct_method.m and in lines 23, 24, 25 enter the feedback gain returned from GA, and Ce, Cu returned from Lemma3.m, and run.
5) The last returns the satisfaction frequencies verified over 10000 disturbance sequence samples (change as needed in line 59) and the solve time of MOSEK solver. It also prints Fig. 1 (left).
