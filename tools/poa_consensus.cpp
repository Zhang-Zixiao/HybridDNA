#include <iostream>
#include <string>
#include "spoa/spoa.hpp"

int main() {
    auto engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 5, -4, -8, -6);
    spoa::Graph graph;
    std::string sequence;
    bool any = false;
    while (std::getline(std::cin, sequence)) {
        if (!sequence.empty() && sequence.back() == '\r') sequence.pop_back();
        if (sequence.empty()) continue;
        auto alignment = engine->Align(sequence, graph);
        graph.AddAlignment(alignment, sequence);
        any = true;
    }
    if (!any) return 1;
    std::cout << graph.GenerateConsensus() << '\n';
    return 0;
}
