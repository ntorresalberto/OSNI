#ifndef TWISTINTEGRATOR_H
#define TWISTINTEGRATOR_H


#include "ode_ab.h"


class TwistIntegrator : public OSNI::ODEAb

{
public:

    TwistIntegrator(const std::shared_ptr<const StrainBasedParameterisation> t_strain_parameterisation) : OSNI::ODEAb(6, t_strain_parameterisation) {}

    TwistIntegrator(const std::shared_ptr<const StrainBasedParameterisation> t_strain_parameterisation,
                         const unsigned int t_number_of_Chebyshev_points) : OSNI::ODEAb(4, t_strain_parameterisation, t_number_of_Chebyshev_points) {}

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point) override;

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point) override;


private:

    Eigen::Vector3d m_xi_at_point { Eigen::Vector3d::Zero() };

    Eigen::Matrix4d m_A_at_point { Eigen::MatrixXd::Zero(6, 6) };

    Eigen::Matrix4d m_b_at_point { Eigen::VectorXd::Zero(6) };
};

#endif // TWISTINTEGRATOR_H
