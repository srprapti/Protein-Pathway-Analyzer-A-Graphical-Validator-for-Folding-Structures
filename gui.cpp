#include <SFML/Graphics.hpp>
#include "gui.hpp"
#include <iostream>
#include <cmath>

using namespace std;

void launchGUI(const vector<char>& sequence,
               const vector<vector<pair<int, int>>>& graph,
               const vector<int>& path,
               const vector<string>& properties,
               const vector<pair<int, int>>& misfoldEdges) {

    int N = sequence.size();
    vector<sf::Vector2f> positions;
    for (int i = 0; i < N; ++i) {
        positions.push_back(sf::Vector2f(100.f + i * 100, 300.f));
    }

    sf::RenderWindow window(sf::VideoMode(1000, 600), "Protein Pathway Visualizer");

    sf::Font font;
    if (!font.loadFromFile("LiberationSans-Regular.ttf")) {
        cout << "Font file (LiberationSans-Regular.ttf) not found." << endl;
        return;
    }

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();

            if (event.type == sf::Event::MouseButtonPressed) {
                sf::Vector2f mouse(event.mouseButton.x, event.mouseButton.y);
                for (int i = 0; i < N; ++i) {
                    float dx = mouse.x - positions[i].x;
                    float dy = mouse.y - positions[i].y;
                    if (sqrt(dx * dx + dy * dy) <= 30.f) {
                        cout << properties[i] << endl;
                    }
                }
            }
        }

        window.clear(sf::Color::White);

        // Draw edges
        for (int u = 0; u < N; ++u) {
            for (auto [v, _] : graph[u]) {
                sf::Color color = sf::Color::Black;

                for (auto [mu, mv] : misfoldEdges) {
                    if (u == mu && v == mv) {
                        color = sf::Color::Red;
                        break;
                    }
                }

                for (size_t i = 1; i < path.size(); ++i) {
                    if (path[i - 1] == u && path[i] == v) {
                        color = sf::Color::Green;
                        break;
                    }
                }

                sf::Vertex line[] = {
                    sf::Vertex(positions[u], color),
                    sf::Vertex(positions[v], color)
                };
                window.draw(line, 2, sf::Lines);
            }
        }

        // Draw nodes
        for (int i = 0; i < N; ++i) {
            sf::CircleShape node(25);
            node.setFillColor(sf::Color(100, 149, 237));
            node.setPosition(positions[i].x - 25, positions[i].y - 25);
            window.draw(node);

            sf::Text label;
            label.setFont(font);
            label.setString(string(1, sequence[i]) + "(" + to_string(i + 1) + ")");
            label.setCharacterSize(14);
            label.setFillColor(sf::Color::White);
            label.setPosition(positions[i].x - 20, positions[i].y - 10);
            window.draw(label);
        }

        window.display();
    }
}
