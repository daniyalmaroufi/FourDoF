string joint_name       # "ITR", "OTR" (rotational_controller) or "ITT", "OTT" (translational_controller)
bool relative            # true = move by `position` (delta), false = move to `position` (absolute)
float64 position         # degrees (rotational) or mm (translational)
float64 velocity         # deg/s or mm/s; 0 = use this joint's current default speed
float64 acceleration     # deg/s^2 or mm/s^2; 0 = use this joint's current default acceleration
