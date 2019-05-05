
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <memory>
#include <random>

#define N_ITERATIONS 100
#define N_PARTICULES 1024

struct QuadTree
{
    protected:
        float m_av_mass;
        float m_x_center;
        float m_y_center;
        bool hasChildren;

        std::vector<QuadTree*> m_children;
        QuadTree();

};

class Particles : QuadTree
{
    private:
        std::vector<float> m_x;
        std::vector<float> m_y;
        std::vector<float> m_vx;
        std::vector<float> m_vy;
        std::vector<float> m_ax;
        std::vector<float> m_ay;
        std::vector<float> m_mass;

    public:
        Particles();
        void initialize(std::string pattern);
        void resetTree();
        void buildTree();
};