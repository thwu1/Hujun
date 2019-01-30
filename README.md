This is Hujun's final project

Main function : vcycle_dgs & vcycle_uzawa

Main change:

GS iterate just one time during a DGS iteration

Rewrite cal_res(the previous one is probably wrong)

What to do next:

Rewrite vcylce, try to restrict the residue(done)

New file:

vcycle_dgs2 & vcycle_uzawa2

performs much better than old ones

UPDATE:

Update vcycle_uzawa2 (The previous one has a mistake)

(vcycle_dgs2 is probably wrong too, need to check later)

NEW FILE:

folder:  explicit_vcycle

RUN demo to TEST (Parameter a can be larger?)

1.30 update:

Run problem1 problem2 problem3 problem4 to test speed
