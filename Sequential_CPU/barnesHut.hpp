
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <memory>
#include <random>
#include <limits>
#include <cmath>

//#define SAVE true

#define THETA 0.5
#define N_ITERATIONS 100
#define N_PARTICULES 65536
#define CHILD 4
#define G 1.0
#define EPSILON 0.0001
#define MASS 1.0

#define BOUNDS 7
#define TIMESTEP 0.001
#define CIRCLE_OFFSET 3
#define FAR_SPACE 30

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
        
        
        std::vector<QuadTree*> m_children;

        QuadTree();
        void quadtreeReset();
        void createNode(int quadrant, double mass, double x, double y, double depth);
        void addBodyToNode(double mass, double x, double y);
        void computeBranchesComponent(double x, double y, double m, double &fx, double &fy);

};

class Particles
{
    private:
        std::vector<double> m_x;
        std::vector<double> m_y;
        std::vector<double> m_vx;
        std::vector<double> m_vy;
        std::vector<double> m_ax;
        std::vector<double> m_ay;
        std::vector<double> m_mass;

        QuadTree m_tree;

        double m_x_min, m_x_max, m_y_min, m_y_max;

        void computePosition(double x, double y, BoxLimits &limits, bool updateLimits);
    public:
        Particles();
        void initialize(std::string pattern);
        void resetTree();
        void buildTree();
        void computeDisplacement();
        #ifdef SAVE
            void saveToFile(std::ofstream *file);
        #endif
};
