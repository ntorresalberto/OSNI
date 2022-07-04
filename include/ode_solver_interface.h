/**
 * @file ode_solver_interface.h
 * @author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * @brief This files contains the base class for integrating a linear ode with the spectral numerical integration
 * @version 0.1
 * @date 03-07-2022
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef ODE_SOLVER_INTERFACE_H
#define ODE_SOLVER_INTERFACE_H

#include <Eigen/Dense>

#include <unsupported/Eigen/KroneckerProduct>
#include "chebyshev_differentiation.hpp"


namespace OSNI {



/*!
 * \brief The SpectralIntegrationBase class
 */
class ODESolverInterface
{
public:
    ODESolverInterface();

    virtual ~ODESolverInterface();


    virtual void integrate(const Eigen::VectorXd& t_initial_state)=0;

    virtual Eigen::VectorXd getStateAtPoint(const unsigned int t_point)const=0;

    inline virtual Eigen::MatrixXd getStack()const=0;

};

}   //  namespace OSNI

#endif // ODE_SOLVER_INTERFACE_H
