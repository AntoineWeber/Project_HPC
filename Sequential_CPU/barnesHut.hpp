
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <memory>
#include <random>
#include <limits>

#define N_ITERATIONS 100
#define N_PARTICULES 2048
#define CHILD 4

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
        
        std::vector<QuadTree*> m_children;
        QuadTree();
        void quadtreeReset();

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
        void createNode(QuadTree *curr_node, int quadrant, float mass, float x, float y);
        void addBodyToNode(QuadTree *curr_node, float mass, float x, float y);

    public:
        Particles();
        void initialize(std::string pattern);
        void resetTree();
        void computeBoundingBox();
        void buildTree();
        void computeForce();
};