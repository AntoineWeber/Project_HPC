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
    int quadrant_internal_node;
    int quadrant_internal_point;

    bool done = false;
    bool constructInternalNode = false;

    BoxLimits limits;
    QuadTree *curr_node;

    for(unsigned int i=0; i<N_PARTICULES; i++)
    {
        // pointer points to root
        curr_node = &m_tree;
        // initialize border of grid
        limits.left = m_x_min;
        limits.right = m_x_max;
        limits.top = m_y_max;
        limits.bottom = m_y_min;

        done = false;
        // descend along the tree until found the right branch
        while(!done)
        {
            // update the current node with the particle information
            addBodyToNode(curr_node, m_mass[i], m_x[i], m_y[i]);

            //locate the particle
            computePosition(m_x[i], m_y[i], limits, true);
            quadrant = limits.quadrant;

            // if no child at this node, create it
            if (curr_node->m_children[quadrant] == nullptr)
            {
                createNode(curr_node, quadrant, m_mass[i], m_x[i], m_y[i]);
                done = true;
            }
            // if already a child, go deeper
            else
            {
                if (curr_node->m_children[quadrant]->hasChildren)
                {
                    // already internal node, go deeper.
                    curr_node = curr_node->m_children[quadrant];
                }
                else
                {
                    // has to bring the body down
                    constructInternalNode = false;
                    while(!constructInternalNode)
                    {
                        computePosition(m_x[i], m_y[i], limits, false);
                        quadrant_internal_point = limits.quadrant;
                        computePosition(curr_node->m_children[quadrant]->m_x_center, \
                                         curr_node->m_children[quadrant]->m_y_center, limits, false);

                        quadrant_internal_node = limits.quadrant;
                        if (quadrant_internal_point != quadrant_internal_node)
                        {
                            createNode(curr_node->m_children[quadrant], quadrant_internal_point, m_mass[i], m_x[i], m_y[i]);
                            createNode(curr_node->m_children[quadrant], quadrant_internal_node, curr_node->m_children[quadrant]->m_av_mass, \
                                         curr_node->m_children[quadrant]->m_x_center, \
                                         curr_node->m_children[quadrant]->m_y_center);
                            addBodyToNode(curr_node->m_children[quadrant], m_mass[i], m_x[i], m_y[i]);

                            constructInternalNode = true;
                            done = true;
                        } 
                        else
                        {
                            computePosition(m_x[i], m_y[i], limits, true);
                            quadrant_internal_point = limits.quadrant;
                            createNode(curr_node->m_children[quadrant], quadrant_internal_point, m_mass[i], m_x[i], m_y[i]);
                            addBodyToNode(curr_node->m_children[quadrant]->m_children[quadrant_internal_point], \
                                            curr_node->m_children[quadrant]->m_av_mass, \
                                            curr_node->m_children[quadrant]->m_x_center,
                                            curr_node->m_children[quadrant]->m_y_center);
                            curr_node = curr_node->m_children[quadrant];
                        } 
                    }
                }
            }
        }
    }
}

void Particles::createNode(QuadTree *curr_node, int quadrant, float mass, float x, float y)
{
    curr_node->m_children[quadrant] = new QuadTree();
    curr_node->m_children[quadrant]->m_av_mass = mass;
    curr_node->m_children[quadrant]->m_x_center = x;
    curr_node->m_children[quadrant]->m_y_center = y;
}

void Particles::addBodyToNode(QuadTree *curr_node, float mass, float x, float y)
{
    curr_node->m_x_center += x*mass / curr_node->m_av_mass;
    curr_node->m_y_center += y*mass / curr_node->m_av_mass;
    curr_node->m_x_center = curr_node->m_x_center * curr_node->m_av_mass / (curr_node->m_av_mass + mass);
    curr_node->m_y_center = curr_node->m_y_center * curr_node->m_av_mass / (curr_node->m_av_mass + mass);
    curr_node->m_av_mass += mass;
    curr_node->hasChildren = true;
}

void Particles::computePosition(float x, float y, BoxLimits &limits, bool updateLimits)
{
    if (x < (limits.left+limits.right)/2)
    {
        if (y >= (limits.top+limits.bottom)/2)
        {
            limits.quadrant = 0;
            if (updateLimits)
            {
                limits.right = (limits.left+limits.right)/2;
                limits.top = (limits.top+limits.bottom)/2;
            }
        }
        else if (y < (limits.top+limits.bottom)/2)
        {
            limits.quadrant = 2;
            if (updateLimits)
            {
                limits.right = (limits.left+limits.right)/2;
                limits.top = (limits.top+limits.bottom)/2;
            }
        }
    }
    else if (x >= (limits.left+limits.right)/2)
    {
        if (y >= (limits.top+limits.bottom)/2)
        {
            limits.quadrant = 1;
            if (updateLimits)
            {
                limits.left = (limits.left+limits.right)/2;
                limits.bottom = (limits.top+limits.bottom)/2;
            }
        }
        else if (y < (limits.top+limits.bottom)/2)
        {
            limits.quadrant = 3;
            if (updateLimits)
            {
                limits.left = (limits.left+limits.right)/2;
                limits.top = (limits.top+limits.bottom)/2;
            }
        }
    }
}