#ifndef GUI_HPP
#define GUI_HPP

#include <vector>
#include <string>

void launchGUI(const std::vector<char>& sequence,
               const std::vector<std::vector<std::pair<int, int>>>& graph,
               const std::vector<int>& shortestPath,
               const std::vector<std::string>& properties,
               const std::vector<std::pair<int, int>>& misfoldEdges);

#endif
