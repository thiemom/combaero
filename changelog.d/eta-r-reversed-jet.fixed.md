- A starved or reversed impingement jet row is now crossflow-controlled:
  CR-3936's `eta_r` is taken at the crossflow over the jet's magnitude, so it
  tends to 1 from either side of zero jet flow (flagged `eta_r_extrapolated`).
  A reversed row used to drop to `eta_r = 0` (the reference temperature fell to
  the jet temperature), unflagged. Against CR-3630 Table D.2 end to end, the
  bypass `eta` RMSE falls from 0.163 to 0.078 (#467).
