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

#include <ur_ik_solver/ur_ik_solver.h>
#include <ur_kinematics/ur_kin.h>
#include <pluginlib/class_list_macros.h>

PLUGINLIB_EXPORT_CLASS(ik_solver::UrIkSolver, ik_solver::IkSolver)


namespace ik_solver
{

inline bool UrIkSolver::customConfig()
{
  if (base_frame_.find("base_link") == std::string::npos)
  {
      ROS_ERROR("%s/base_frame should be set equal to [PREFIX]base_link instead of %s",nh_.getNamespace().c_str(),base_frame_.c_str());
      return false;
  }
  if (flange_frame_.find("tool0") == std::string::npos)
  {
      ROS_ERROR("%s/flange_frame should be set equal to [PREFIX]tool0 instead of %s",nh_.getNamespace().c_str(),flange_frame_.c_str());
      return false;
  }

  Eigen::AngleAxisd link6_ee(0.5*M_PI,Eigen::Vector3d::UnitZ());
  Eigen::AngleAxisd link6_tool0(-0.5*M_PI,Eigen::Vector3d::UnitX());


  T_flange_ee_.setIdentity();
  T_flange_ee_.linear()=link6_tool0.toRotationMatrix().inverse()*link6_ee.toRotationMatrix();
  return true;
}


std::vector<Eigen::VectorXd> UrIkSolver::getIk(const Eigen::Affine3d& T_base_flange,
                                                   const std::vector<Eigen::VectorXd> & seeds,
                                                   const int& desired_solutions,
                                                   const int& max_stall_iterations)
{
  double q_sols_array[n_sol*n_joints];
  std::vector<Eigen::VectorXd> q_sols;
  // From Affine3d Column-major to row-major

  Eigen::Affine3d T_base_ee=T_base_flange*T_flange_ee_;
  Eigen::Transform<double,3,Eigen::Affine,Eigen::RowMajor> T_base_flange_rm = T_base_ee;

  int q_sols_found = ur_kinematics::inverse(T_base_flange_rm.data(), q_sols_array);

  for(int idx = 0; idx < q_sols_found; idx++)
  {
    Eigen::VectorXd single_sol = Eigen::Map<Eigen::VectorXd>(q_sols_array+idx*n_joints,n_joints,1);
    std::vector<Eigen::VectorXd> multiple_sols = getMultiplicity(single_sol);
    for(const Eigen::VectorXd& tmp : multiple_sols)
    {
      if(outOfBound(tmp))
        continue;
      q_sols.push_back(tmp);
    }
  }

  return q_sols;
}

inline std::vector<Eigen::VectorXd> UrIkSolver::getMultiplicity(const Eigen::VectorXd &q)
{
  std::vector<std::vector<double>> multiturn_ax(n_joints);

  for (unsigned int idx = 0; idx < n_joints; idx++)
  {
    multiturn_ax.at(idx).push_back(q(idx));
    double tmp=q(idx);
    while (true)
    {
      tmp+=2*M_PI;
      if (tmp>ub_(idx))
        break;
      multiturn_ax.at(idx).push_back(tmp);
    }
    tmp=q(idx);
    while (true)
    {
      tmp-=2*M_PI;
      if (tmp<lb_(idx))
        break;
      multiturn_ax.at(idx).push_back(tmp);
    }
  }

  std::vector<Eigen::VectorXd> multiturn;
  multiturn.push_back(q);
  for (unsigned int idx = 0; idx < n_joints; idx++)
  {
    size_t size_multiturn=multiturn.size();
    for (size_t is=1;is<multiturn_ax.at(idx).size();is++)
    {
      for (size_t im=0;im<size_multiturn;im++)
      {
        Eigen::VectorXd new_q=multiturn.at(im);
        new_q(idx)=multiturn_ax.at(idx).at(is);
        multiturn.push_back(new_q);
      }
    }
  }

  return multiturn;

}

}   // end namespace ik_solver
