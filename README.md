# ur_ik_solver

Analytic inverse kinematics plugin for `ik_solver` (UR3, UR3e, UR5, UR5e, UR10, UR10e, UR16e, UR20,
UR30).

## URDF refinement

The closed-form solver uses the nominal DH parameters. A URDF generated with the robot calibration
(`kinematics_parameters_file` of `ur_description`) differs from them, and the analytic solutions can be
some mm away from the requested pose (up to about 3 mm on a UR10e).

At configuration the plugin builds the chain `base_frame` -> `flange_frame` from the URDF. Each analytic
solution is then refined with Gauss-Newton iterations on that chain (SVD pseudo-inverse of the Jacobian,
step limited to 0.2 rad, position and rotation error below 1e-10), and `getFK` uses the same chain, so
IK and FK agree with the URDF. Solutions that do not converge are discarded, and two that converge to
the same configuration are kept once.

Near a singularity two nominal solutions merge, while the URDF can have two distinct ones about 0.1 rad
apart, or the pose can be slightly outside the nominal workspace. Extra seeds are used there:

- arguments of acos/asin out of [-1, 1] by less than 0.05 are clamped to the border;
- if s_min/s_max of the Jacobian is below 0.05 (shoulder, elbow or wrist singularity), seeds at
  0.05 and 0.15 rad on both sides along the null direction (right singular vector of s_min);
- if |sin(q5)| < 0.2 (wrist), 7 more analytic solutions with q6 = 2*pi*k/8, along the nominal self-motion.

On 100 poses sampled near each singularity (UR10e with the Sharework calibration), the configuration that
generated the pose is found in 96-100 % of the cases; with the plain analytic seeds and a damped step it
was 20-35 %. Away from singularities only the 8 analytic seeds are refined.

| parameter | default | meaning |
|---|---|---|
| `refine_with_urdf` | `true` | refine on the URDF chain; `false`: nominal DH kinematics only |

If the chain cannot be built (frames not in the URDF, a joint that is not revolute, a joint of the
chain missing in `joint_names`), the plugin logs a warning and uses the nominal kinematics only.
