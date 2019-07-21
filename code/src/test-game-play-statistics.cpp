//
// Created by Mikael Zayenz Lagerkvist
//


#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <chrono>
#include <random>

#include <gecode/driver.hh>
#include <gecode/int.hh>
#include <gecode/minimodel.hh>

#include "patchwork/lib.h"

const double avg(const std::vector<int>& v) {
    double sum = 0;
    for (const auto &item : v) {
        sum += item;
    }
    return sum / v.size();
}



int main(int argc, char **argv) {
    const int iterations = 100;

    std::vector<int> green_plies;
    std::vector<int> yellow_plies;
    std::vector<int> alternatives;
    std::vector<int> alternatives_branching;
    std::vector<int> score_diff;
    double green_wins = 0;
    double yellow_wins = 0;
    double draws = 0;

    for (int i = 0; i < iterations; ++i) {
        std::cerr << "Playing game " << i << " of " << iterations << std::endl;
        int green_moves = 0;
        int yellow_moves = 0;
        patchwork::Context context = patchwork::Context(42 + i)
                .with(patchwork::Mode::MODE_SINGLE,
                      patchwork::PlacementPolicy::PP_BOTTOM_LEFT_ALL_ROTATIONS,
                      patchwork::PlacementEvaluation::PE_MINIMUM_PLACEMENT_REGRET);
        patchwork::Context alt_context = context
                .with(patchwork::Mode::MODE_BRANCHING,
                        patchwork::PlacementPolicy::PP_ALL,
                        patchwork::PlacementEvaluation::PE_FIRST);
        patchwork::State *current = patchwork::State::create(context);
        while (!current->finished()) {
            //current->print(std::cerr);
            if (current->current_player() == patchwork::Player::GREEN) {
                ++green_moves;
            } else {
                ++yellow_moves;
            }
            {
                int alts = 0;
                std::vector<patchwork::State *> branching_test = current->next_states(alt_context, alts);
                alternatives_branching.emplace_back(alts);
                for (auto & l : branching_test) {
                    delete l;
                }
            }

            int alts = 0;
            std::vector<patchwork::State *> next_states = current->next_states(context, alts);
            alternatives.emplace_back(alts);
            assert(next_states.size() == 1);

            delete current;
            current = next_states[0];
        }


        int green_score = current->current_score(patchwork::Player::GREEN);
        int yellow_score = current->current_score(patchwork::Player::YELLOW);
        if (green_score > yellow_score) {
            ++green_wins;
        } else if (green_score < yellow_score) {
            ++yellow_wins;
        } else {
            ++draws;
        }
        score_diff.emplace_back(abs(green_score - yellow_score));

        delete current;
        green_plies.emplace_back(green_moves);
        yellow_plies.emplace_back(yellow_moves);
    }

    std::cout << "Green wins " << green_wins << " yellow wins " << yellow_wins << " draws " << draws << std::endl;
    std::cout << "Average absolute score diff " << avg(score_diff) << std::endl;
    std::cout << "Average Green moves " << avg(green_plies) << std::endl;
    std::cout << "Average Yellow moves " << avg(yellow_plies) << std::endl;
    std::cout << "Average alterantives " << avg(alternatives) << std::endl;
    std::cout << "Average branchings " << avg(alternatives_branching) << std::endl;

    return 0;
}
