#include "barnesHut.hpp"


/**
Constructor of the class Quadtree. Initialize every value to 0 and the pointers to the 4 quadrants to
null.
*/
Node::Node()
{
    m_av_mass = 0.0;
    m_x_center = 0.0;
    m_y_center = 0.0;
    m_s = 0.0;
    hasChildren = false;
}

/**
Constructor for the class particle.
*/
Particles::Particles(int n_nodes)
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

    m_tree = thrust::host_vector<Node>(n_nodes, Node());
}

/**
To call the quadtree reset from the particle class.
*/
void Particles::resetTree(int n_nodes)
{
    m_tree.clear();
    m_tree = thrust::host_vector<Node>(n_nodes, Node());
}

/**
Method building the barnes-hut tree
*/
void Particles::buildTree()
{   // The quadrant variable can take value from 0 to 3
    // 0:NW   1:NE  2:SW    3:SE 

    int quadrant = -1;
    int quadrant_internal_point = -1;
    int quadrant_internal_node = -1;

    bool done = false;
    bool constructInternalNode = false;

    double end_leaf_mass,end_leaf_posx,end_leaf_posy;

    int absoluteOffset = 0;
    int depthOffset = 0;
    int nDepthNode = 1;
    int depth = 0;

    // reset tree after having deleted everything
    m_tree[0].m_s = (double) 2*BOUNDS;
    m_tree[0].hasChildren = true;

    BoxLimits limits;

    for(unsigned int i=0; i<N_PARTICULES; i++)
    {
        // if particle in the far space
        if (m_x[i] > FAR_SPACE || m_x[i] < -FAR_SPACE || m_y[i] > FAR_SPACE || m_y[i] < -FAR_SPACE)
        {
            continue;
        }

        // initialize border of grid
        limits.left = m_x_min;
        limits.right = m_x_max;
        limits.top = m_y_max;
        limits.bottom = m_y_min;

        // add particle to root node
        if (i == 0)
        {
            createNode(0, m_x[i], m_y[i], m_mass[i]);
        }
        else
        {
            addBodyToNode(0, m_x[i], m_y[i], m_mass[i]);
        }
        
        depth = 1;
        done = false;
        depthOffset = 0;
        absoluteOffset = 1;
        nDepthNode = 4;
        // descend along the tree until found the right branch
        while(!done)
        {
            if (depth == MAX_DEPTH-1)
            {
                break;
            }

            //locate the particle
            computePosition(m_x[i], m_y[i], limits, true);
            quadrant = limits.quadrant;
            
            // if first layer, move depthOffset to the quadrant to be on the right node
            if (absoluteOffset == 1)
            {
                depthOffset = quadrant;
            }
            // else add the quadrant to be on the right node
            else
            {
                depthOffset += quadrant;
            }
            
            // if no child at this node, create it and go to next particle
            if (m_tree[absoluteOffset + depthOffset].m_av_mass == 0)
            {
                createNode(absoluteOffset + depthOffset, m_x[i], m_y[i], m_mass[i]);
                done = true;
            }
            // if already exists, go deeper down the tree
            else
            {
                if (m_tree[absoluteOffset + depthOffset].hasChildren == true)
                {
                    // add the point to the node
                    addBodyToNode(absoluteOffset + depthOffset, m_x[i], m_y[i], m_mass[i]);

                    // update to the start of the children. Specific quadrant will be updated beginning
                    // of the next while iteration.
                    updateOffsets(absoluteOffset, depthOffset, nDepthNode, depth);
                }
                else
                {
                    // node has no children, hence is not an internal node. Have to create the rest of the branch and bring both
                    // body down : the particle considered and the considered node being a particle.
                    constructInternalNode = false;
                    end_leaf_mass = m_tree[absoluteOffset + depthOffset].m_av_mass;
                    end_leaf_posx = m_tree[absoluteOffset + depthOffset].m_x_center;
                    end_leaf_posy = m_tree[absoluteOffset + depthOffset].m_y_center;
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
                        if (depth == MAX_DEPTH-1)
                        {
                            // add particule before quitting
                            addBodyToNode(absoluteOffset + depthOffset + quadrant_internal_node,  m_x[i], m_y[i], m_mass[i]);
                            break;
                        }
                        computePosition(m_x[i], m_y[i], limits, false);
                        quadrant_internal_point = limits.quadrant;
                        computePosition(end_leaf_posx, end_leaf_posy, limits, false);
                        quadrant_internal_node = limits.quadrant;
                        // if after having cut the space they both go to a different quadrant, construct both nodes and go to next particle
                        if (quadrant_internal_point != quadrant_internal_node)
                        {
                            // add the point to the parent node
                            addBodyToNode(absoluteOffset+depthOffset, m_x[i], m_y[i], m_mass[i]);
                            m_tree[absoluteOffset+depthOffset].hasChildren = true;

                            // create the two children
                            updateOffsets(absoluteOffset, depthOffset, nDepthNode, depth);

                            createNode(absoluteOffset + depthOffset + quadrant_internal_point, m_x[i], m_y[i], m_mass[i]);
                            createNode(absoluteOffset + depthOffset + quadrant_internal_node, end_leaf_posx, end_leaf_posy, end_leaf_mass);

                            
                            constructInternalNode = true;
                            done = true;
                        } 
                        else
                        {
                            // add the point to the parent node
                            addBodyToNode(absoluteOffset+depthOffset, m_x[i], m_y[i], m_mass[i]);
                            m_tree[absoluteOffset+depthOffset].hasChildren = true;

                            // else create internal nodes until they go to different quadrants.
                            // need to call this function to cut the space in 4
                            computePosition(end_leaf_posx, end_leaf_posy, limits, true);
                            quadrant_internal_node = limits.quadrant;

                            // create children in the same state as the entering of the while loop
                            // (end_leaf being at the end) and iterate until you can separate them
                            // or went to max depth.
                            updateOffsets(absoluteOffset, depthOffset, nDepthNode, depth);
                            depthOffset += quadrant_internal_node;

                            createNode(absoluteOffset + depthOffset, end_leaf_posx, end_leaf_posy, end_leaf_mass);
                        } 
                    }
                }
            }
        }
    }
}

void Particles::createNode(int offset, double x, double y, double m)
{
    m_tree[offset].m_x_center = x;
    m_tree[offset].m_y_center = y;
    m_tree[offset].m_av_mass = m;
    m_tree[offset].hasChildren = false;
}

void Particles::addBodyToNode(int offset, double x, double y, double m)
{
    m_tree[offset].m_x_center += x * m / m_tree[offset].m_av_mass;
    m_tree[offset].m_y_center += y * m / m_tree[offset].m_av_mass;

    m_tree[offset].m_x_center = m_tree[offset].m_x_center * m_tree[offset].m_av_mass / (m_tree[offset].m_av_mass + m);
    m_tree[offset].m_y_center = m_tree[offset].m_y_center * m_tree[offset].m_av_mass / (m_tree[offset].m_av_mass + m);

    m_tree[offset].m_av_mass += m;
}

void updateOffsets(int &absolOff, int &depthOff, int &nNode, int &depth)
{
    depthOff = depthOff * 4;
    absolOff += nNode;
    nNode *= 4;
    depth += 1;
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
        std::cout << "x is " << x << " and y is " << y <<std::endl;
        std::cout << "limits are : up : " << limits.top << " down : " << limits.bottom << std::endl;
        std::cout << " left : " << limits.left << " right : " << limits.right << std::endl;
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