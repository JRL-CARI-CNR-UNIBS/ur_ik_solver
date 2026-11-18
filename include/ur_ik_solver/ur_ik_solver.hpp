/*
Copyright (c) 2022, JRL-CARI CNR-STIIMA/UNIBS
Manuel Beschi manuel.beschi@unibs.it
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef UR_IK_SOLVER__UR_IK_SOLVER
#define UR_IK_SOLVER__UR_IK_SOLVER

#include <ik_solver/ik_solver.hpp>
#include <Eigen/Geometry>

namespace ik_solver
{
class UrIkSolver: public IkSolver
{
public:
  virtual Solutions getIk(const Eigen::Affine3d& T_base_flange,
                                     const std::vector<Eigen::VectorXd> & seeds,
                                     const int& desired_solutions,
                                     const int& min_stall_iterations = -1,
                                     const int& max_stall_iterations = -1) override;

  virtual Eigen::Affine3d getFK(const Configuration& s) override;
  virtual bool config(const ParamMap& params) override;
protected:
  /** One URDF joint between base_frame and flange_frame. */
  struct ChainJoint
  {
    Eigen::Affine3d origin;  ///< Pose of the joint frame in the parent link.
    Eigen::Vector3d axis;    ///< Rotation axis in the joint frame.
    int index;               ///< Position in joint_names_, -1 for a fixed joint.
  };

  /** Builds chain_ from the URDF; false if the chain is missing or has non-revolute joints. */
  bool buildChain();
  /** Forward kinematics of the URDF chain (base_frame -> flange_frame). */
  Eigen::Affine3d chainFK(const Configuration& q) const;
  /** Forward kinematics and geometric Jacobian (base frame) of the URDF chain. */
  Eigen::Affine3d chainJacobian(const Configuration& q, Eigen::Matrix<double, 6, 6>& J) const;
  /**
   * Newton iterations on the URDF chain, starting from an analytic solution. The analytic solver
   * uses the nominal DH parameters, the URDF can contain the robot calibration (errors of some mm).
   * Near a singularity the result can be far from the analytic solution.
   * Returns false (q unchanged) if it does not converge.
   */
  bool refine(Configuration& q, const Eigen::Affine3d& T_base_flange) const;
  /**
   * Seeds of refine(): the analytic solutions (arguments slightly out of the acos/asin domain clamped)
   * plus, near a singularity, configurations along the null direction of the Jacobian (smallest
   * singular value) and, near the wrist singularity, along the nominal wrist self-motion.
   */
  std::vector<Configuration> refinementSeeds(const double* T_base_ee_row_major) const;

  Eigen::Affine3d T_flange_ee_;
  Eigen::Affine3d T_ee_flange_;
  const unsigned int n_joints = 6;
  const unsigned int n_sol = 8;
  std::vector<ChainJoint> chain_;  ///< Empty: analytic kinematics only.
};
}  //  namespace ik_solver

#endif // UR_IK_SOLVER__UR_IK_SOLVER
