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

#include "ur_ik_solver/ur_ik_solver.hpp"
#include "ur_ik_solver/ur_kinematics/ur_kin.h"

#include <algorithm>
#include <array>
#include <Eigen/SVD>
#include <urdf_model/model.h>

#if ROS_X == 1
#include <pluginlib/class_list_macros.h>
#else
#include <pluginlib/class_list_macros.hpp>
#endif

PLUGINLIB_EXPORT_CLASS(ik_solver::UrIkSolver, ik_solver::IkSolver)

namespace ik_solver
{

bool UrIkSolver::config(const ParamMap& params)
{
  if(!IkSolver::config(params))
  {
    IK_FATAL(this->logger_, "Ik solver initial config FAILED");
    return false;
  }
  IK_INFO(this->logger_, "Check joint number");
  if(this->joint_names_.size() != 6)
  {
    IK_FATAL(this->logger_, "UR Ik Solver found " << this->joint_names_.size() << " joints instead of 6");
    return false;
  }
  Eigen::AngleAxisd link6_ee(0.5*M_PI,Eigen::Vector3d::UnitZ());
  Eigen::AngleAxisd link6_tool0(-0.5*M_PI,Eigen::Vector3d::UnitX());

  T_flange_ee_.setIdentity();
  T_flange_ee_.linear()=link6_tool0.toRotationMatrix().inverse()*link6_ee.toRotationMatrix();
  T_ee_flange_ = T_flange_ee_.inverse();
  bool refine_with_urdf = true;
  get_optional_param(this->logger_, params_, "refine_with_urdf", refine_with_urdf);
  chain_.clear();
  if (!refine_with_urdf)
    IK_INFO(this->logger_, "UR Ik Solver: refine_with_urdf is false, nominal DH kinematics only");
  else if (buildChain())
    IK_INFO(this->logger_, "UR Ik Solver: analytic solutions refined on the URDF chain " << base_frame_ << " -> "
                                                                                          << flange_frame_);
  else
    IK_WARN(this->logger_, "UR Ik Solver: cannot use the URDF chain, nominal DH kinematics only");
  IK_INFO(this->logger_, "UR Ik Solver ready");
  return true;
}

bool UrIkSolver::buildChain()
{
  chain_.clear();
  if (!model_)
    return false;
  std::vector<ChainJoint> reversed;
  std::vector<bool> found(joint_names_.size(), false);
  urdf::LinkConstSharedPtr link = model_->getLink(flange_frame_);
  while (link && link->name != base_frame_)
  {
    const urdf::JointSharedPtr& joint = link->parent_joint;
    if (!joint)
      return false;  // base_frame is not an ancestor of flange_frame
    const urdf::Pose& o = joint->parent_to_joint_origin_transform;
    ChainJoint c;
    c.origin = Eigen::Translation3d(o.position.x, o.position.y, o.position.z) *
               Eigen::Quaterniond(o.rotation.w, o.rotation.x, o.rotation.y, o.rotation.z);
    c.axis = Eigen::Vector3d(joint->axis.x, joint->axis.y, joint->axis.z);
    c.index = -1;
    if (joint->type == urdf::Joint::REVOLUTE || joint->type == urdf::Joint::CONTINUOUS)
    {
      auto it = std::find(joint_names_.begin(), joint_names_.end(), joint->name);
      if (it == joint_names_.end())
        return false;
      c.index = static_cast<int>(it - joint_names_.begin());
      found[c.index] = true;
    }
    else if (joint->type != urdf::Joint::FIXED)
      return false;
    reversed.push_back(c);
    link = model_->getLink(joint->parent_link_name);
  }
  if (!link || std::find(found.begin(), found.end(), false) != found.end())
    return false;
  chain_.assign(reversed.rbegin(), reversed.rend());
  return true;
}

Eigen::Affine3d UrIkSolver::chainFK(const Configuration& q) const
{
  Eigen::Affine3d T = Eigen::Affine3d::Identity();
  for (const ChainJoint& c : chain_)
  {
    T = T * c.origin;
    if (c.index >= 0)
      T = T * Eigen::AngleAxisd(q(c.index), c.axis);
  }
  return T;
}

Eigen::Affine3d UrIkSolver::chainJacobian(const Configuration& q, Eigen::Matrix<double, 6, 6>& J) const
{
  // Geometric Jacobian in the base frame
  Eigen::Affine3d T = Eigen::Affine3d::Identity();
  std::vector<std::pair<int, std::pair<Eigen::Vector3d, Eigen::Vector3d>>> axes;  // index, (axis, point)
  for (const ChainJoint& c : chain_)
  {
    T = T * c.origin;
    if (c.index >= 0)
    {
      axes.push_back({ c.index, { T.linear() * c.axis, T.translation() } });
      T = T * Eigen::AngleAxisd(q(c.index), c.axis);
    }
  }
  J.setZero();
  for (const auto& [index, axis_point] : axes)
  {
    const auto& [z, p] = axis_point;
    J.block<3, 1>(0, index) += z.cross(T.translation() - p);
    J.block<3, 1>(3, index) += z;
  }
  return T;
}

std::vector<Configuration> UrIkSolver::refinementSeeds(const double* T_base_ee_row_major) const
{
  // Calibration errors move a pose by some mm: near the border of the nominal workspace the
  // arguments of acos/asin leave [-1, 1] by about 1e-2
  constexpr double domain_tolerance = 0.05;
  // Near a shoulder, elbow or wrist singularity two analytic solutions merge in the nominal model,
  // while the URDF can have two distinct ones, apart by about 0.1 rad along the null direction
  constexpr double singular_ratio = 0.05;  // s_min / s_max of the Jacobian
  constexpr std::array<double, 2> null_offsets = { 0.05, 0.15 };  // [rad]
  // Near q5 = 0 or pi the nominal solutions form a self-motion of the wrist, longer than the offsets
  constexpr double near_wrist = 0.2;  // |sin(q5)|
  constexpr int wrist_samples = 8;    // q6 values along the nominal self-motion

  double q_array[n_sol * n_joints];
  std::vector<Configuration> seeds;
  bool near_wrist_singularity = false;
  const int found = ur_kinematics::inverse(T_base_ee_row_major, q_array, 0.0, domain_tolerance);
  for (int idx = 0; idx < found; idx++)
  {
    Configuration q = Eigen::Map<Eigen::VectorXd>(q_array + idx * n_joints, n_joints, 1);
    if (!q.allFinite())
      continue;
    seeds.push_back(q);
    near_wrist_singularity |= std::abs(std::sin(q(4))) < near_wrist;

    Eigen::Matrix<double, 6, 6> J;
    chainJacobian(q, J);
    Eigen::JacobiSVD<Eigen::Matrix<double, 6, 6>> svd(J, Eigen::ComputeFullV);
    const Eigen::Matrix<double, 6, 1>& s = svd.singularValues();
    if (s(5) < singular_ratio * s(0))
    {
      const Eigen::Matrix<double, 6, 1> v = svd.matrixV().col(5);
      for (double offset : null_offsets)
        for (double sign : { -1.0, 1.0 })
          seeds.push_back(q + sign * offset * v);
    }
  }
  if (near_wrist_singularity)
  {
    for (int k = 1; k < wrist_samples; k++)
    {
      const int n = ur_kinematics::inverse(T_base_ee_row_major, q_array, 2.0 * M_PI * k / wrist_samples,
                                           domain_tolerance, near_wrist);
      for (int idx = 0; idx < n; idx++)
      {
        Configuration q = Eigen::Map<Eigen::VectorXd>(q_array + idx * n_joints, n_joints, 1);
        if (q.allFinite())
          seeds.push_back(q);
      }
    }
  }
  return seeds;
}

bool UrIkSolver::refine(Configuration& q, const Eigen::Affine3d& T_base_flange) const
{
  constexpr int max_iterations = 50;
  constexpr double tolerance = 1e-10;  // [m] and [rad]
  constexpr double max_step = 0.2;     // [rad] the seeds near a singularity can be far from the solution
  constexpr double rcond = 1e-12;      // singular values below rcond * s_max are dropped

  // Gauss-Newton with the SVD pseudo-inverse. A damped least squares step would be cut along the
  // direction of the smallest singular value, the one that matters near a singularity.
  Configuration x = q;
  Eigen::Matrix<double, 6, 6> J;
  for (int iter = 0; iter < max_iterations; ++iter)
  {
    const Eigen::Affine3d T = chainJacobian(x, J);
    Eigen::Matrix<double, 6, 1> e;
    e.head<3>() = T_base_flange.translation() - T.translation();
    Eigen::AngleAxisd rot(T_base_flange.linear() * T.linear().transpose());
    e.tail<3>() = rot.angle() * rot.axis();
    if (e.head<3>().norm() < tolerance && e.tail<3>().norm() < tolerance)
    {
      q = x;
      return true;
    }
    Eigen::JacobiSVD<Eigen::Matrix<double, 6, 6>> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Eigen::Matrix<double, 6, 1>& s = svd.singularValues();
    Eigen::Matrix<double, 6, 1> projected = svd.matrixU().transpose() * e;
    for (int i = 0; i < 6; i++)
      projected(i) = s(i) > rcond * s(0) ? projected(i) / s(i) : 0.0;
    Configuration dx = svd.matrixV() * projected;
    const double largest = dx.cwiseAbs().maxCoeff();
    if (largest > max_step)
      dx *= max_step / largest;
    x += dx;
  }
  return false;
}

Solutions UrIkSolver::getIk(const Eigen::Affine3d& T_base_flange,
                                                   const Configurations& seeds,
                                                   const int& desired_solutions,
                                                   const int& /*min_stall_iterations*/,
                                                   const int& /*max_stall_iterations*/)
{
  double q_sols_array[n_sol*n_joints];
  Configurations q_sols;
  // From Affine3d Column-major to row-major

  Eigen::Affine3d T_base_ee=T_base_flange*T_flange_ee_;
  Eigen::Transform<double,3,Eigen::Affine,Eigen::RowMajor> T_base_flange_rm = T_base_ee;

  std::vector<Eigen::VectorXd> tmp_sols;
  if (chain_.empty())
  {
    int q_sols_found = ur_kinematics::inverse(T_base_flange_rm.data(), q_sols_array);
    for (int idx = 0; idx < q_sols_found; idx++)
      tmp_sols.push_back(Eigen::Map<Eigen::VectorXd>(q_sols_array + idx * n_joints, n_joints, 1));
  }
  else
  {
    for (Configuration q : refinementSeeds(T_base_flange_rm.data()))
    {
      if (!refine(q, T_base_flange))
        continue;  // not a solution of the URDF kinematics
      // Seeds near a singularity can converge to the same configuration
      const bool duplicate = std::any_of(tmp_sols.begin(), tmp_sols.end(), [&q](const Configuration& other) {
        const Eigen::ArrayXd d = (other - q).array() / (2.0 * M_PI);
        return ((d - d.round()) * 2.0 * M_PI).abs().maxCoeff() < 1e-5;
      });
      if (!duplicate)
        tmp_sols.push_back(q);
    }
  }
  q_sols = getMultiplicity(tmp_sols, ub(), lb(), {1,1,1,1,1,1});

  for (auto it = q_sols.begin(); it != q_sols.end(); )
  {
      std::vector<int> oob = outOfBound(*it, ub(), lb());
      // erase any solution that has at least one out-of-bound joint
      if (!oob.empty())
          it = q_sols.erase(it);
      else
          ++it;
  }

  Solutions sols;
  sols.configurations() = q_sols;
  sols.message() = "Found " + std::to_string(q_sols.size()) + " solutions";
  return sols;
}

Eigen::Affine3d UrIkSolver::getFK(const Configuration& s)
{
  if (!chain_.empty())
    return chainFK(s);
  double T[16];
  ur_kinematics::forward(s.data(), T); // fills a row-major 4x4

  // Map as a row-major matrix type (note: RowMajor is part of the Matrix type, not Map Options)
  Eigen::Map<const Eigen::Matrix<double,4,4,Eigen::RowMajor>> H_rm(T);

  Eigen::Affine3d T_base_ee(H_rm);           // Eigen handles storage; semantics: t in last column

  Eigen::Affine3d T_base_flange=T_base_ee*T_ee_flange_;

  return T_base_flange;
}

}   // end namespace ik_solver
