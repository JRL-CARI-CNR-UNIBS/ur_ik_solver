^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Changelog for package ur_ik_solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Forthcoming
-----------
* Refine the analytic solutions on the URDF chain (damped Newton) and compute getFK on it, so that a
  calibrated URDF is matched exactly (the nominal DH parameters were up to about 3 mm off on a UR10e).
  Parameter refine_with_urdf (default true) restores the previous behaviour.
* Near singularities: acos/asin arguments slightly out of domain are clamped, extra seeds along the
  Jacobian null direction and along the nominal wrist self-motion, SVD pseudo-inverse step.
* ur_kinematics::inverse: optional arguments domain_tolerance and wrist_singular_tolerance (defaults keep
  the previous behaviour).

0.1.0 (2026-10-22)
------------------
* Fix build type forced to Debug: CMAKE_BUILD_TYPE was hard-coded in CMakeLists.txt, overriding
  -DCMAKE_BUILD_TYPE. The code was optimized only when cnr_logger was installed as Release (it exports
  -Ofast to its consumers), otherwise it compiled at -O0. The package now defaults to Release
  (-O3 -DNDEBUG) and honours -DCMAKE_BUILD_TYPE. Rebuild with colcon build --cmake-clean-cache.
* Contributors: Federico Parma, Manuel Beschi, Nicola Pedrocchi, Roberto Fausti
