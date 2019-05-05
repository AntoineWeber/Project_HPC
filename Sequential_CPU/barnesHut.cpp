#include "barnesHut.hpp"


QuadTree::QuadTree()
{
    m_av_mass = 0;
    m_x_center = 0;
    m_y_center = 0;
    hasChildren = false;

Particles::Particles()
{
    m_x.reserve(N_PARTICULES);
    m_y.reserve(N_PARTICULES);
    m_vx.reserve(N_PARTICULES);
    m_vy.reserve(N_PARTICULES);
    m_ax.reserve(N_PARTICULES);
    m_ay.reserve(N_PARTICULES);
    m_mass.reserve(N_PARTICULES);
}

void Particles::initialize(std::string pattern)
{
    if (pattern == "uniform")
    {
        std::cout << "Initializing the particles in a square with uniform distribution" << std::endl;
        std::default_random_engine generator(time(0));
        std::uniform_real_distribution<float> square(-5.0, 5.0);

        // loop through all particles
        for (int i = 0; i < N_PARTICULES; i++){

            m_mass[i] = 1.0;
            m_x[i] = square(generator);
            m_y[i] = square(generator);

            // set velocity and acceleration at 0
            m_vx[i] = 0;
            m_vy[i] = 0;
            m_ax[i] = 0;
            m_ay[i] = 0;
        }
    }
}

void Particles::resetTree()
{
    m_children.clear();
    for (unsigned int i=0; i<m_children.size(); i++)
    {
        m_children[i]->Reset();
    }
}

void Particles::buildTree()
{

}