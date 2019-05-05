#include "barnesHut.hpp"


quadTree::quadTree()
{
    // instanciate quadtree for children
    m_children = std::make_unique<quadTree[]>(4);
}