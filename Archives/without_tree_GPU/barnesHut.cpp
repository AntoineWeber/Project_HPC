#include "barnesHut.hpp"


/**
Constructor of the class Quadtree. Initialize every value to 0 and the pointers to the 4 quadrants to
null.
*/
QuadTree::QuadTree()
{
    m_av_mass = 0.0;
    m_x_center = 0.0;
    m_y_center = 0.0;
    m_s = 0.0;
    hasChildren = false;
    m_children.resize(CHILD);
    for (unsigned int i=0; i<CHILD; i++)
    {
        m_children[i] = nullptr;
    }
}

/**
Method used to reset the tree. Recursive calls to all children
*/
void QuadTree::quadtreeReset()
{
    for (unsigned int i=0; i<CHILD; i++)
    {
        if (m_children[i] != nullptr)
        {
            m_children[i]->quadtreeReset();
        }
        delete m_children[i];
    }
    m_children.clear();
    hasChildren = false;
}

d_Node::d_Node()
{
    m_av_mass = 0.0;
    m_x_center = 0.0;
    m_y_center = 0.0;
    m_s = 0.0;
    hasChildren = false;

    m_children = new d_Node*[CHILD];
    for (unsigned int i=0; i<CHILD; i++)
    {
        m_children[i] = nullptr;
    }
}

/**
Constructor for the class particle.
*/
Particles::Particles()
{
    m_x.resize(N_PARTICULES);
    m_y.resize(N_PARTICULES);
    m_vx.resize(N_PARTICULES);
    m_vy.resize(N_PARTICULES);
    m_ax.resize(N_PARTICULES);
    m_ay.resize(N_PARTICULES);
    m_mass.resize(N_PARTICULES);

    m_x_min = -FAR_SPACE;
    m_x_max = FAR_SPACE;
    m_y_min = -FAR_SPACE;
    m_y_max = FAR_SPACE;
}

/**
To call the quadtree reset from the particle class.
*/
void Particles::resetTree()
{
    m_tree.quadtreeReset();
}


/*
void Particles::computeBoundingBox()
{
    for (unsigned int i=0; i<N_PARTICULES; i++)
    {
        if (m_x[i] > m_x_max)
        {
            m_x_max = m_x[i];
        }
        else if (m_x[i] < m_x_min)
        {
            m_x_min = m_x[i];
        }
        
        if (m_y[i] > m_y_max)
        {
            m_y_max = m_y[i];
        }
        else if (m_y[i] < m_y_min)
        {
            m_y_min = m_y[i];
        }
    }
}
*/

/**
Method building the barnes-hut tree
*/
void Particles::buildTree()
{   // The quadrant variable can take value from 0 to 3
    // 0:NW   1:NE  2:SW    3:SE 

    int quadrant;
    int quadrant_internal_node;
    int quadrant_internal_point;

    bool done = false;
    bool constructInternalNode = false;

    BoxLimits limits;
    limits.quadrant = -1;
    QuadTree *curr_node;

    double end_leaf_mass, end_leaf_posx, end_leaf_posy;

    // reset tree after having deleted everything
    m_tree.m_s = (double) 2*BOUNDS;
    m_tree.m_children.resize(CHILD);
    for (unsigned int i=0; i<CHILD; i++)
    {
        m_tree.m_children[i] = nullptr;
    }

    for(unsigned int i=0; i<N_PARTICULES; i++)
    {
        // if particle in the far space
        if (m_x[i] > FAR_SPACE || m_x[i] < -FAR_SPACE || m_y[i] > FAR_SPACE || m_y[i] < -FAR_SPACE)
        {
            continue;
        }
        // pointer points to root of the tree
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
            curr_node->addBodyToNode(m_mass[i], m_x[i], m_y[i]);

            //locate the particle
            computePosition(m_x[i], m_y[i], limits, true);
            quadrant = limits.quadrant;

            // if no child at this node, create it and go to next particle
            if (curr_node->m_children[quadrant] == nullptr)
            {
                curr_node->createNode(quadrant, m_mass[i], m_x[i], m_y[i], curr_node->m_s);
                done = true;
            }
            // if already a child, go deeper down the tree
            else
            {
                if (curr_node->m_children[quadrant]->hasChildren)
                {
                    // node has already children, go deeper
                    curr_node = curr_node->m_children[quadrant];
                }
                else
                {
                    // node has no children, hence is not an internal node. Have to create the rest of the branch and bring both
                    // body down : the particle considered and the considered node being a particle.
                    constructInternalNode = false;
                    end_leaf_mass = curr_node->m_children[quadrant]->m_av_mass;
                    end_leaf_posx = curr_node->m_children[quadrant]->m_x_center;
                    end_leaf_posy = curr_node->m_children[quadrant]->m_y_center;

                    // if tree not destroyed properly, without the forcing computation,
                    // can get stuck in an infinite loop here trying to separate the same particle.
                    // if really not lucky, the randomization may have created twice the same particle.
                    if (end_leaf_posx == m_x[i] && end_leaf_posy == m_y[i])
                    {
                        std::cout << "Same particle ! Try re-running." << std::endl;
                        exit(0);
                    }

                    while(!constructInternalNode)
                    {
                        computePosition(m_x[i], m_y[i], limits, false);
                        quadrant_internal_point = limits.quadrant;
                        computePosition(end_leaf_posx, end_leaf_posy, limits, false);
                        quadrant_internal_node = limits.quadrant;
                        // if after having cut the space they both go to a different quadrant, construct both nodes and go to next particle
                        if (quadrant_internal_point != quadrant_internal_node)
                        {
                            curr_node->m_children[quadrant]->createNode(quadrant_internal_point, m_mass[i], m_x[i], m_y[i],
                                                                         curr_node->m_children[quadrant]->m_s);
                            curr_node->m_children[quadrant]->createNode(quadrant_internal_node, end_leaf_mass, end_leaf_posx,
                                                                         end_leaf_posy, curr_node->m_children[quadrant]->m_s);
                            curr_node->m_children[quadrant]->addBodyToNode(m_mass[i], m_x[i], m_y[i]);

                            constructInternalNode = true;
                            done = true;
                        } 
                        else
                        {
                            curr_node->m_children[quadrant]->addBodyToNode(m_mass[i], m_x[i], m_y[i]);
                            // else create internal nodes until they go to different quadrants.
                            // need to call this function to cut the space in 4
                            computePosition(end_leaf_posx, end_leaf_posy, limits, true);
                            quadrant_internal_point = limits.quadrant;
                            curr_node->m_children[quadrant]->createNode(quadrant_internal_point, end_leaf_mass, end_leaf_posx, end_leaf_posy,
                                                                         curr_node->m_children[quadrant]->m_s);
                            curr_node = curr_node->m_children[quadrant];
                            quadrant = quadrant_internal_point;
                        } 
                    }
                }
            }
        }
    }
}

/**
Method used to create a node when none exists in the tree
*/
void QuadTree::createNode(int quadrant, double mass, double x, double y, double prof)
{
    m_children[quadrant] = new QuadTree();
    m_children[quadrant]->m_av_mass = mass;
    m_children[quadrant]->m_x_center = x;
    m_children[quadrant]->m_y_center = y;
    m_children[quadrant]->m_s = prof / 2.0;
    //std::cout << m_children[quadrant]->m_s << std::endl;
}

/**
Method used to add a body to an already existing node
*/
void QuadTree::addBodyToNode(double mass, double x, double y)
{
    m_x_center += x*mass / m_av_mass;
    m_y_center += y*mass / m_av_mass;
    m_x_center = m_x_center * m_av_mass / (m_av_mass + mass);
    m_y_center = m_y_center * m_av_mass / (m_av_mass + mass);
    m_av_mass += mass;
    hasChildren = true;
}


/**
Method used to locate the particle in the space defined by limits, and update the borders if boolean to true.
As the function takes a reference as input, there is no need to return any value.
*/
void Particles::computePosition(double x, double y, BoxLimits &limits, bool updateLimits)
{
    // if on the left side of the space
    if (x < (limits.left+limits.right)/2)
    {
        // if on the top side
        if (y >= (limits.top+limits.bottom)/2)
        {
            // NW quadrant
            limits.quadrant = 0;
            if (updateLimits)
            {
                limits.right = (limits.left+limits.right)/2;
                limits.top = (limits.top+limits.bottom)/2;
            }
        }
        // if on the bottom side
        else if (y < (limits.top+limits.bottom)/2)
        {
            // SW quadrant
            limits.quadrant = 2;
            if (updateLimits)
            {
                limits.right = (limits.left+limits.right)/2;
                limits.top = (limits.top+limits.bottom)/2;
            }
        }
    }
    // if on the right side of the space
    else if (x >= (limits.left+limits.right)/2)
    {
        // if on the top side
        if (y >= (limits.top+limits.bottom)/2)
        {
            // quadrant NE
            limits.quadrant = 1;
            if (updateLimits)
            {
                limits.left = (limits.left+limits.right)/2;
                limits.bottom = (limits.top+limits.bottom)/2;
            }
        }
        // if on the bottom side
        else if (y < (limits.top+limits.bottom)/2)
        {
            // quadrant SE
            limits.quadrant = 3;
            if (updateLimits)
            {
                limits.left = (limits.left+limits.right)/2;
                limits.top = (limits.top+limits.bottom)/2;
            }
        }
    }
    else
    {
        std::cout << "Should not happen" << std::endl;
        exit(0);
    }
    
}

#ifdef SAVE
void Particles::saveToFile(std::ofstream *file)
{
    for (unsigned int i=0; i<N_PARTICULES; i++)
    {
        *file << m_x[i] << " " << m_y[i] << "\n";
    }
}
#endif