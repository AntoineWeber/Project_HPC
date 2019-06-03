
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <memory>
#include <random>
#include <limits>
#include <cmath>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

//#define SAVE true

#define THETA 0.5
#define N_ITERATIONS 100
#define CHILD 4
#define G 1.0
#define EPSILON 0.0001
#define MASS 1.0

#define BOUNDS 7
#define TIMESTEP 0.001
#define CIRCLE_OFFSET 3
#define FAR_SPACE 30

#define N_PARTICULES 65536
#define GRID_SIZE 32
#define BLOCK_SIZE 32

struct BoxLimits
{
    double left,right,top,bottom;
    int quadrant;
};

class QuadTree
{
    public:
        double m_av_mass;
        double m_x_center;
        double m_y_center;
        double m_s;

        bool hasChildren;
        
        thrust::host_vector<QuadTree*> m_children;

        QuadTree();
        void quadtreeReset();
        void createNode(int quadrant, double mass, double x, double y, double depth);
        void addBodyToNode(double mass, double x, double y);
};

class Particles
{
    public:
        thrust::host_vector<double> m_x;
        thrust::host_vector<double> m_y;
        thrust::host_vector<double> m_vx;
        thrust::host_vector<double> m_vy;
        thrust::host_vector<double> m_ax;
        thrust::host_vector<double> m_ay;
        thrust::host_vector<double> m_mass;

        QuadTree m_tree;

        double m_x_min, m_x_max, m_y_min, m_y_max;

        void computePosition(double x, double y, BoxLimits &limits, bool updateLimits);

        Particles();
        void initialize(std::string pattern);
        void resetTree();
        void buildTree();
        void computeDisplacement();
        #ifdef SAVE
            void saveToFile(std::ofstream *file);
        #endif
};

class d_Node
{
    public:
        double m_av_mass;
        double m_x_center;
        double m_y_center;
        double m_s;

        bool hasChildren;
        
        d_Node** m_children;

        d_Node();
};
