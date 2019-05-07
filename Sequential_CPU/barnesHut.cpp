#include "barnesHut.hpp"


QuadTree::QuadTree()
{
    m_av_mass = 0;
    m_x_center = 0;
    m_y_center = 0;
    hasChildren = false;
    m_children.reserve(CHILD);
    for (unsigned int i=0; i<CHILD; i++)
    {
        m_children[i] = nullptr;
    }
}

void QuadTree::quadtreeReset()
{
    for (unsigned int i=0; i<m_children.size(); i++)
    {
        m_children[i]->quadtreeReset();
        delete m_children[i];
    }
    m_children.clear();
    hasChildren = false;
}

Particles::Particles()
{
    m_x.reserve(N_PARTICULES);
    m_y.reserve(N_PARTICULES);
    m_vx.reserve(N_PARTICULES);
    m_vy.reserve(N_PARTICULES);
    m_ax.reserve(N_PARTICULES);
    m_ay.reserve(N_PARTICULES);
    m_mass.reserve(N_PARTICULES);

    m_x_min = std::numeric_limits<float>::infinity();
    m_x_max = -std::numeric_limits<float>::infinity();
    m_y_min = std::numeric_limits<float>::infinity();
    m_y_max = -std::numeric_limits<float>::infinity();
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
    m_tree.quadtreeReset();
}

// might get rid of this function if considering a fixed box and everything outside is "far space"
void Particles::computeBoundingBox()
{
    for (unsigned int i=0; i<N_PARTICULES; i++)
    {
        if (m_x[i] < m_x_min)
        {
            m_x_min = m_x[i];
        }
        else if (m_x[i] > m_x_max)
        {
            m_x_max = m_x[i];
        }
        
        if (m_y[i] < m_y_min)
        {
            m_y_min = m_y[i];
        }
        else if (m_y[i] > m_y_max)
        {
            m_y_max = m_y[i];
        }
    }
}

void Particles::buildTree()
{
    // 0:NW   1:NE  2:SW    3:SE  
    int quadrant;
    bool done = false;
    float l = m_x_min;
    float r = m_x_max;
    float t = m_y_max;
    float d = m_y_min;

    QuadTree *curr_node = &m_tree;

    for(unsigned int i=0; i<N_PARTICULES; i++)
    {
        // first update the current node with the particle information
        curr_node->m_x_center += m_x[i]*m_mass[i] / curr_node->m_av_mass;
        curr_node->m_y_center += m_y[i]*m_mass[i] / m_tree.m_av_mass;
        curr_node->m_x_center = curr_node->m_x_center * curr_node->m_av_mass / (curr_node->m_av_mass + m_mass[i]);
        curr_node->m_y_center = curr_node->m_y_center * curr_node->m_av_mass / (curr_node->m_av_mass + m_mass[i]);
        curr_node->m_av_mass += m_mass[i];
        curr_node->hasChildren = true;

        done = false;
        // descend along the tree until found the right branch
        while(!done)
        {
            //locate the particle
            if (m_x[i] < (l+r)/2)
            {
                if (m_y[i] > (t+d)/2)
                {
                    quadrant = 0; 
                    r = (l+r)/2;
                    d = (t+d)/2;
                }
                else if (m_y[i] < (t+d)/2)
                {
                    quadrant = 2;
                    r = (l+r)/2;
                    t = (t+d)/2;
                }
            }
            else if (m_x[i] > (l+r)/2)
            {
                if (m_y[i] > (t+d)/2)
                {
                    quadrant = 1;
                    l = (l+r)/2;
                    d = (t+d)/2;
                }
                else if (m_y[i] < (t+d)/2)
                {
                    quadrant = 3;
                    l = (l+r)/2;
                    t = (t+d)/2;
                }

            }

            // if no child at this node, update it
            if (curr_node->m_children[quadrant] == nullptr)
            {
                curr_node->m_children[quadrant] = new QuadTree();
                curr_node->m_children[quadrant]->m_av_mass = m_mass[i];
                curr_node->m_children[quadrant]->m_x_center = m_x[i];
                curr_node->m_children[quadrant]->m_y_center = m_y[i];
                done = true;
            }

            else
            {
                // if leaf already exists, update it
                curr_node = curr_node->m_children[quadrant];
                /*
                curr_node->m_x_center += m_x[i]*m_mass[i] / curr_node->m_av_mass;
                curr_node->m_y_center += m_y[i]*m_mass[i] / m_tree.m_av_mass;
                curr_node->m_x_center = curr_node->m_x_center * curr_node->m_av_mass / (curr_node->m_av_mass + m_mass[i]);
                curr_node->m_y_center = curr_node->m_y_center * curr_node->m_av_mass / (curr_node->m_av_mass + m_mass[i]);
                curr_node->m_av_mass += m_mass[i];
                curr_node->hasChildren = true;
                */
            }
            
        }

    }
}