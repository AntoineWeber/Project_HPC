#include "barnesHut.hpp"



/**
Constructor of the class Quadtree. Initialize every value to 0
*/
__host__ Node::Node()
{
    m_av_mass = 0.0;
    m_x_center = 0.0;
    m_y_center = 0.0;
    m_s = 0.0;
    hasChildren = false;

    children = 0;
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

    m_tree = std::make_unique<Node[]>(n_nodes);
}

/**
Method to clean the tree and allocate a new one
*/
void Particles::resetTree(int n_nodes)
{
    m_tree.reset();
    m_tree = std::make_unique<Node[]>(n_nodes);
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

    double prev_prof;

    BoxLimits limits;

    // loop on all particles
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
            createNode(0,0,1, m_x[i], m_y[i], m_mass[i], (double) 2*2*FAR_SPACE);
            m_tree[0].hasChildren = true;
        }
        else
        {
            addBodyToNode(0, m_x[i], m_y[i], m_mass[i]);
        }

        prev_prof = m_tree[0].m_s;
        depth = 1;
        done = false;
        // absolute offset is the offset to the first particle for a given depth
        // depth offset is the offset the the beginning of the four quadrants, from the absolute offset
        // finally quadrant is the quadrant, being an integer between 0 and 3
        // nDeptNode is the number of nodes for the given depth
        depthOffset = 0;
        absoluteOffset = 1;
        nDepthNode = 4;
        // descend along the tree until found the right branch
        while(!done)
        {
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
                createNode(absoluteOffset, depthOffset, nDepthNode, m_x[i], m_y[i], m_mass[i], prev_prof);
                done = true;
            }
            // if already exists, go deeper down the tree
            else
            {
                // if child exists and already have children
                if (m_tree[absoluteOffset + depthOffset].hasChildren == true)
                {
                    // add the point to the node
                    addBodyToNode(absoluteOffset + depthOffset, m_x[i], m_y[i], m_mass[i]);

                    // update to the start of the children. Specific quadrant will be updated beginning
                    // of the next while iteration.
                    updateOffsets(absoluteOffset, depthOffset, nDepthNode, depth, prev_prof);
                }
                else
                {
                    // node has no children, hence is not an internal node. Have to create the rest of the branch and bring both
                    // body down : the particle considered and the considered node being a particle.
                    constructInternalNode = false;
                    end_leaf_mass = m_tree[absoluteOffset + depthOffset].m_av_mass;
                    end_leaf_posx = m_tree[absoluteOffset + depthOffset].m_x_center;
                    end_leaf_posy = m_tree[absoluteOffset + depthOffset].m_y_center;

                    // if really not lucky, the randomization may have created twice the same particle.
                    if (end_leaf_posx == m_x[i] && end_leaf_posy == m_y[i])
                    {
                        std::cout << "Same particle ! Try re-running." << std::endl;
                        exit(0);
                    }

                    while(!constructInternalNode)
                    {
                        // if arrived at the max depth, add particle to last layer and leave
                        if (depth == MAX_DEPTH-1)
                        {
                            // add particule before quitting
                            addBodyToNode(absoluteOffset + depthOffset,  m_x[i], m_y[i], m_mass[i]);
                            done = true;
                            break;
                        }
                        // locate both the external node and the current particle
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
                            updateOffsets(absoluteOffset, depthOffset, nDepthNode, depth, prev_prof);

                            createNode(absoluteOffset, depthOffset + quadrant_internal_point, nDepthNode, m_x[i], m_y[i], m_mass[i], prev_prof);
                            createNode(absoluteOffset, depthOffset + quadrant_internal_node, nDepthNode, end_leaf_posx, end_leaf_posy, end_leaf_mass, prev_prof);

                            
                            constructInternalNode = true;
                            done = true;
                        } 
                        else
                        {
                            // If both particle and external node go deeper in the same quadrant, iterate until they 
                            // go to different quadrants.

                            // add the point to the parent node
                            addBodyToNode(absoluteOffset+depthOffset, m_x[i], m_y[i], m_mass[i]);
                            m_tree[absoluteOffset+depthOffset].hasChildren = true;

                            // need to call this function to cut the space in 4
                            computePosition(end_leaf_posx, end_leaf_posy, limits, true);
                            quadrant_internal_node = limits.quadrant;

                            // create children in the same state as the entering of the while loop
                            // (end_leaf being at the end) and iterate until you can separate them
                            // or went to max depth.
                            updateOffsets(absoluteOffset, depthOffset, nDepthNode, depth, prev_prof);
                            depthOffset += quadrant_internal_node;

                            createNode(absoluteOffset, depthOffset, nDepthNode, end_leaf_posx, end_leaf_posy, end_leaf_mass, prev_prof);
                        } 
                    }
                }
            }
        }
    }
}

/**
 * Method to create a node when no particle is already on it
 */
void Particles::createNode(int absOff, int depthOff, int nNode, double x, double y, double m, double prof)
{
    m_tree[absOff + depthOff].m_x_center = x;
    m_tree[absOff + depthOff].m_y_center = y;
    m_tree[absOff + depthOff].m_av_mass = m;
    m_tree[absOff + depthOff].hasChildren = false;
    // m_s represents the size of the current quadtree. To be used when computing the
    // ratio s/r to determine if the node is far enough in the compute force part.
    m_tree[absOff + depthOff].m_s = prof / 2.0;

    // offset to the 4 children of the node
    m_tree[absOff + depthOff].children = absOff+nNode + depthOff*4;
}
/**
 * Method to add a particle on an already existing node.
 * Updates the center of mass and the total mass
 */
void Particles::addBodyToNode(int offset, double x, double y, double m)
{
    m_tree[offset].m_x_center += x * m / m_tree[offset].m_av_mass;
    m_tree[offset].m_y_center += y * m / m_tree[offset].m_av_mass;

    m_tree[offset].m_x_center = m_tree[offset].m_x_center * m_tree[offset].m_av_mass / (m_tree[offset].m_av_mass + m);
    m_tree[offset].m_y_center = m_tree[offset].m_y_center * m_tree[offset].m_av_mass / (m_tree[offset].m_av_mass + m);

    m_tree[offset].m_av_mass += m;
}

/**
 * Function updating the offsets. The purpose of this function is to bring the offsets
 * to the beginning of the children of the considered node.
 * It also updates the depth, the number of nodes on this depth and the size of the space of the quadtrees for this depth
 */
void updateOffsets(int &absolOff, int &depthOff, int &nNode, int &depth, double &prev_prof)
{
    depthOff = depthOff * 4;
    absolOff += nNode;
    nNode *= 4;
    depth += 1;
    prev_prof = prev_prof / 2.0;
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
        // If particle not within the space hence in the far space. This
        // method should never be called for such particles.
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