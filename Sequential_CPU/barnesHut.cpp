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
Method used to initialize the particles.
*/
void Particles::initialize(std::string pattern)
{
    if (pattern == "uniform")
    {
        std::cout << "Initializing the particles in a square with uniform distribution" << std::endl;
        std::default_random_engine generator(time(0));
        std::uniform_real_distribution<double> square(-BOUNDS, BOUNDS);

        // loop through all particles
        for (int i = 0; i < N_PARTICULES; i++){

            m_mass[i] = MASS;
            m_x[i] = square(generator);
            m_y[i] = square(generator);

            // set velocity and acceleration at 0
            m_vx[i] = 0;
            m_vy[i] = 0;
            m_ax[i] = 0;
            m_ay[i] = 0;
        }
    }
    else if (pattern == "circle")
    {
        std::cout << "Initializing the particles in a square with a circle distribution" << std::endl;
        std::default_random_engine generator;
        std::uniform_real_distribution<double> radius(CIRCLE_OFFSET, 2*BOUNDS);
        std::uniform_real_distribution<double> angle(0, 2*M_PI);
        double r, theta;
        // loop through all particles
        for (int i = 0; i < N_PARTICULES; i++){
            /*
            if (i==0)
            {
                m_mass[i] = MASS * 10000;
                m_x[i] = 0.0;
                m_y[i] = 0.0;
                m_vx[i] = 0.0;
                m_vy[i] = 0.0;
                m_ax[i] = 0.0;
                m_ay[i] = 0.0;

            }
            else
            {
            */
            r = radius(generator);
            theta = angle(generator);

            m_mass[i] = MASS;
            m_x[i] = r*cos(theta);
            m_y[i] = r*sin(theta);

            // set velocity and acceleration at 0
            m_vx[i] = 0;//0.006*sin(theta);
            m_vy[i] = 0;//-0.006*cos(theta);
            m_ax[i] = 0;
            m_ay[i] = 0;
            //}
        }
    }
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

    // reset root after having cleaned the tree.
    m_tree.m_s = (double) 2*FAR_SPACE;
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

                    // if really not lucky, the randomization may have created twice the same particle.
                    if (end_leaf_posx == m_x[i] && end_leaf_posy == m_y[i])
                    {
                        std::cout << "Same particle ! Try re-running." << std::endl;
                        exit(0);
                    }

                    while(!constructInternalNode)
                    {
                        // compute the quadrants of the considered particle and external node
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
                            // else create internal nodes until they go to different quadrants.
                            curr_node->m_children[quadrant]->addBodyToNode(m_mass[i], m_x[i], m_y[i]);

                            // recall this function to cut the space in 4
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
Method used to create a node when none exists 
*/
void QuadTree::createNode(int quadrant, double mass, double x, double y, double prof)
{
    m_children[quadrant] = new QuadTree();
    m_children[quadrant]->m_av_mass = mass;
    m_children[quadrant]->m_x_center = x;
    m_children[quadrant]->m_y_center = y;
    m_children[quadrant]->m_s = prof / 2.0;
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


void Particles::computeDisplacement()
{
    double fx = 0;
    double fy = 0;

    for(unsigned int i=0; i<N_PARTICULES; i++)
    {
        fx = 0;
        fy = 0;

        // if particle in the far space
        if (m_x[i] > FAR_SPACE || m_x[i] < -FAR_SPACE || m_y[i] > FAR_SPACE || m_y[i] < -FAR_SPACE)
        {
            fx = 0;
            fy = 0;
        }
        else
        {
            // call recursive method to compute the force
            m_tree.computeBranchesComponent(m_x[i], m_y[i], m_mass[i], fx, fy);
        }

        // compute acceleration velocity and position based on the force
        m_ax[i] = fx / m_mass[i];
        m_ay[i] = fy / m_mass[i];

        m_vx[i] += m_ax[i]*TIMESTEP;
        m_vy[i] += m_ay[i]*TIMESTEP;

        m_x[i] += 0.5*m_ax[i]*TIMESTEP*TIMESTEP + m_vx[i]*TIMESTEP;
        m_y[i] += 0.5*m_ay[i]*TIMESTEP*TIMESTEP + m_vy[i]*TIMESTEP;
    }
}

void QuadTree::computeBranchesComponent(double x, double y, double m, double &fx, double &fy)
{
    double r, s;
    for (unsigned int i=0; i<m_children.size(); i++)
    {
        // if child exists
        if (m_children[i] != nullptr)
        {
            r = sqrt((m_children[i]->m_x_center-x)*(m_children[i]->m_x_center-x) + (m_children[i]->m_y_center-y)*(m_children[i]->m_y_center-y));
            // if external node, compute force with it and add it to force component
            if (!m_children[i]->hasChildren)
            {
                if (r > EPSILON)
                {
                    fx += (G*m*m_children[i]->m_av_mass*(m_children[i]->m_x_center-x))/r;
                    fy += (G*m*m_children[i]->m_av_mass*(m_children[i]->m_y_center-y))/r;
                }
                else
                {
                    fx += (G*m*m_children[i]->m_av_mass*(m_children[i]->m_x_center-x))/EPSILON;
                    fy += (G*m*m_children[i]->m_av_mass*(m_children[i]->m_y_center-y))/EPSILON;
                }
                
            }
            // if internal node
            else
            {
                // compute the quotient s/r
                s = m_children[i]->m_s;

                // if too far, consider as a single body
                if (s / r < THETA)
                {
                    // still checks for EPSILON as theta could be as small as wanted.
                    if (r > EPSILON)
                    {
                        fx += (G*m*m_children[i]->m_av_mass*(m_children[i]->m_x_center-x))/r;
                        fy += (G*m*m_children[i]->m_av_mass*(m_children[i]->m_y_center-y))/r;
                    }
                    else
                    {
                        fx += (G*m*m_children[i]->m_av_mass*(m_children[i]->m_x_center-x))/EPSILON;
                        fy += (G*m*m_children[i]->m_av_mass*(m_children[i]->m_y_center-y))/EPSILON;
                    }
                }
                else
                {
                    // if close enough, recursive call
                    m_children[i]->computeBranchesComponent(x,y,m,fx,fy);
                }
            }
        }
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