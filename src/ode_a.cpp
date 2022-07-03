#include "ode_a.h"

#include <iostream>


namespace OSNI {


void ODEA::integrate(const Eigen::VectorXd& t_initial_state)
{
    std::cout << "D_IN : \n" << m_D_IN << "\n \n" << std::endl;
    std::cout << "D_NN : \n" << m_D_NN << "\n \n" << std::endl;

    updateA();

    m_ivp = m_D_IN*t_initial_state;
    std::cout << "m_ivp : \n" << m_ivp << "\n \n" << std::endl;

    std::cout << "res : \n" << - m_A_NN.inverse() * m_ivp << "\n \n" << std::endl;

    m_states_stack = - m_A_NN.inverse() * m_ivp;

}




void ODEA::updateA()
{
    for(unsigned int current_chebyshev_point=0; current_chebyshev_point<m_number_of_Chebyshev_points-1; current_chebyshev_point++){



        //  Compute the A matrix of Q' = 1/2 A(K) Q
        m_A_at_Chebyshev_point = computeMatrixAtChebyshevPoint(current_chebyshev_point);

        std::cout << "m_A_at_Chebyshev_point : \n" << m_A_at_Chebyshev_point << "\n \n" << std::endl;

        for (unsigned int row = 0; row < m_state_dimension; ++row) {
            for (unsigned int col = 0; col < m_state_dimension; ++col) {
                int row_index = row*(m_number_of_Chebyshev_points-1) + current_chebyshev_point;
                int col_index = col*(m_number_of_Chebyshev_points-1) + current_chebyshev_point;
                m_A_NN(row_index, col_index) = m_D_NN(row_index, col_index) - m_A_at_Chebyshev_point(row, col);
            }
        }

    }

    std::cout << "m_A_NN : \n" << m_A_NN << "\n \n" << std::endl;


}


}
