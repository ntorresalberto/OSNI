#include "spectral_integration_base.h"


namespace OSNI {

ODESolverInterface::ODESolverInterface(const unsigned int t_state_dimension,
                                                 const std::shared_ptr<const StrainBasedParameterisation> t_strain_parameterisation)
                        : m_state_dimension(t_state_dimension),
                          m_strain_parameterisation(t_strain_parameterisation)  {}


ODESolverInterface::ODESolverInterface(const unsigned int t_state_dimension,
                                                 const std::shared_ptr<const StrainBasedParameterisation> t_strain_parameterisation,
                                                 const unsigned int t_number_of_Chebyshev_points)
                        : m_state_dimension(t_state_dimension),
                          m_strain_parameterisation(t_strain_parameterisation),
                          m_number_of_Chebyshev_points(t_number_of_Chebyshev_points){}

}   //  namespace OSNI
