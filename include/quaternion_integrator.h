#ifndef QUATERNIONINTEGRATOR_H
#define QUATERNIONINTEGRATOR_H

#include "ode_a.h"


class QuaternionIntegrator : public OSNI::ODEA
{
public:
    QuaternionIntegrator(const std::shared_ptr<const StrainBasedParameterisation> t_strain_parameterisation) : OSNI::ODEA(4, t_strain_parameterisation) {}

    QuaternionIntegrator(const std::shared_ptr<const StrainBasedParameterisation> t_strain_parameterisation,
                         const unsigned int t_number_of_Chebyshev_points) : OSNI::ODEA(4, t_strain_parameterisation, t_number_of_Chebyshev_points) {}

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) override;


private:

    Eigen::Vector3d m_K_at_point { Eigen::Vector3d::Zero() };

    Eigen::Matrix4d m_A_at_point { Eigen::Matrix4d::Zero() };

};

#endif // QUATERNIONINTEGRATOR_H
