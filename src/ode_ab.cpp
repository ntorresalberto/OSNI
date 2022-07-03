#include "ode_ab.h"

#include <iostream>


namespace OSNI {


void ODEAb::integrate(const Eigen::VectorXd& t_initial_state)
{

    updateA();
    updateb();

    m_ivp = m_D_IN*t_initial_state;


    m_states_stack = m_A_NN.inverse() * (m_b_NN - m_ivp);

}




void ODEAb::updateA()
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

void ODEAb::updateb()
{
    for(unsigned int current_chebyshev_point=0; current_chebyshev_point<m_number_of_Chebyshev_points-1; current_chebyshev_point++){



        //  Compute the A matrix of Q' = 1/2 A(K) Q
        m_b_at_Chebyshev_point = computerParametersVectorAtPoint(current_chebyshev_point);


        for (unsigned int row = 0; row < m_state_dimension; ++row) {
            m_b_NN(current_chebyshev_point*m_state_dimension + row) = m_b_at_Chebyshev_point(row);
        }

    }


}


}
