
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <memory>


class quadTree
{
    public:
        float m_mass;
        float m_x_center;
        float m_y_center;

        // allocate quadtree for children
        std::unique_ptr<quadTree[]> m_children;

        quadTree();
};