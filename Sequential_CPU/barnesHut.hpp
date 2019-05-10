
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <memory>
#include <random>
#include <limits>
#include <cmath>

#define SAVE true

#define THETA 0.5
#define N_ITERATIONS 100
#define N_PARTICULES 1024
#define CHILD 4
#define G 6.67408e-11
#define EPSILON 0.001
#define MASS 100.0

#define BOUNDS 5
#define TIMESTEP 100

struct BoxLimits
{
    float left,right,top,bottom;
    int quadrant;
};

class QuadTree
{
    public:
        float m_av_mass;
        float m_x_center;
        float m_y_center;
        bool hasChildren;
        int depth;
        
        std::vector<QuadTree*> m_children;

        QuadTree();
        void quadtreeReset();
        void createNode(int quadrant, float mass, float x, float y, int depth);
        void addBodyToNode(float mass, float x, float y);
        void computeBranchesComponent(float x, float y, float m, float &fx, float &fy);

};

class Particles
{
    private:
        std::vector<float> m_x;
        std::vector<float> m_y;
        std::vector<float> m_vx;
        std::vector<float> m_vy;
        std::vector<float> m_ax;
        std::vector<float> m_ay;
        std::vector<float> m_mass;

        QuadTree m_tree;

        float m_x_min, m_x_max, m_y_min, m_y_max;

        void computePosition(float x, float y, BoxLimits &limits, bool updateLimits);
    public:
        Particles();
        void initialize(std::string pattern);
        void resetTree();
        void computeBoundingBox();
        void buildTree();
        void computeDisplacement();
        #ifdef SAVE
            void saveToFile(std::ofstream *file);
        #endif
};