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

class EvaluationStatistics {
private:
    patchwork::PlacementPolicy policy_;
    patchwork::PlacementEvaluation evaluation_;
    int tests_;
    std::vector<double> runtime_;
    double total_runtime_;
    std::vector<int> first_fail_;
    int total_first_fail_;
    std::vector<int> placed_area_first_fail_;
    int total_placed_area_first_fail_;
    std::vector<int> placed_area_;
    int total_placed_area_;
    std::vector<int> placed_tiles_;
    int total_placed_tiles_;
    std::vector<int> failures_;
    int total_failures_;
    std::vector<int> alternatives_;
    int total_alternatives_;
public:
    EvaluationStatistics(const patchwork::PlacementPolicy policy, const patchwork::PlacementEvaluation evaluation)
            : evaluation_(evaluation), policy_(policy),
              tests_(0),
              runtime_{}, total_runtime_(0),
              first_fail_{}, total_first_fail_(0),
              placed_area_first_fail_{}, total_placed_area_first_fail_(0),
              placed_area_{}, total_placed_area_(0),
              placed_tiles_{}, total_placed_tiles_(0),
              failures_{}, total_failures_(0),
              alternatives_{}, total_alternatives_(0)
              {}

    EvaluationStatistics() = default;
    EvaluationStatistics(const EvaluationStatistics& es) = default;

    void add(double runtime, int first_fail, int placed_area_first_fail, int placed_area, int placed_tiles, int failures, const std::vector<int>& alternatives) {
        ++tests_;
        runtime_.emplace_back(runtime);
        total_runtime_ += runtime;
        first_fail_.emplace_back(first_fail);
        total_first_fail_ += first_fail;
        placed_area_first_fail_.emplace_back(placed_area_first_fail);
        total_placed_area_first_fail_ += placed_area_first_fail;
        placed_area_.emplace_back(placed_area);
        total_placed_area_ += placed_area;
        placed_tiles_.emplace_back(placed_tiles);
        total_placed_tiles_ += placed_tiles;
        failures_.emplace_back(failures);
        total_failures_ += failures;
        for (const auto &alternative : alternatives) {
            alternatives_.emplace_back(alternative);
            total_alternatives_ += alternative;
        }
    }

    int tests() const {
        return tests_;
    }

    double mean_runtime() const {
        return total_runtime_ / tests_;
    }

    double mean_runtime_per_alt() const {
        return total_runtime_ / (tests_ * mean_alternatives());
    }

    double mean_first_fail() const {
        return static_cast<double>(total_first_fail_) / tests_;
    }

    double mean_placed_area_first_fail() const {
        return static_cast<double>(total_placed_area_first_fail_) / tests_;
    }

    double mean_placed_area() const {
        return static_cast<double>(total_placed_area_) / tests_;
    }

    int median_placed_area() const {
        if (placed_area_.size() == 1) {
            return placed_area_[0];
        }
        std::vector tmp(placed_area_);
        std::sort(tmp.begin(), tmp.end());
        int median_index = tmp.size() / 2;
        return tmp[median_index];
    }

    int min_placed_area() const {
        return *std::min_element(placed_area_.begin(), placed_area_.end());
    }

    int max_placed_area() const {
        return *std::max_element(placed_area_.begin(), placed_area_.end());
    }

    double mean_placed_tiles() const {
        return static_cast<double>(total_placed_tiles_) / tests_;
    }

    int median_placed_tiles() const {
        if (placed_tiles_.size() == 1) {
            return placed_tiles_[0];
        }
        std::vector tmp(placed_tiles_);
        std::sort(tmp.begin(), tmp.end());
        int median_index = tmp.size() / 2;
        return tmp[median_index];
    }

    int min_placed_tiles() const {
        return *std::min_element(placed_tiles_.begin(), placed_tiles_.end());
    }

    int max_placed_tiles() const {
        return *std::max_element(placed_tiles_.begin(), placed_tiles_.end());
    }

    double mean_failures() const {
        return static_cast<double>(total_failures_) / tests_;
    }

    int max_failures() const {
        return *std::max_element(failures_.begin(), failures_.end());
    }

    double mean_alternatives() const {
        return static_cast<double>(total_alternatives_) / alternatives_.size();
    }

    int max_alternatives() const {
        return *std::max_element(alternatives_.begin(), alternatives_.end());
    }


    std::string to_string() const {
        unsigned long longest_policy = 0;
        for (const auto &policy : patchwork::all_policies) {
            longest_policy = std::max(longest_policy, patchwork::to_string(policy).size());
        }
        unsigned long longest_evaluation = 0;
        for (const auto &evaluation : patchwork::all_evaluations) {
            longest_evaluation = std::max(longest_evaluation, patchwork::to_string(evaluation).size());
        }
        //<< std::right << std::setw(2)
        std::stringstream res;
        res.setf(std::ios::floatfield,std::ios::fixed);
        res << "EvalStat(<"
            << std::right << std::setw(longest_policy) << patchwork::to_string(policy_)
            << ","
            << std::right << std::setw(longest_evaluation) << patchwork::to_string(evaluation_)
            << ">, tests=" << tests_
            << ", avg-runtime=" << std::right << std::setw(9) << std::setprecision(4) << mean_runtime()
            << ", avg-runtime-per-alt=" << std::right << std::setw(9) << std::setprecision(4) << mean_runtime_per_alt()
            << ", avg-ff=" << std::right << std::setw(5) << std::setprecision(2) << mean_first_fail()
            << ", avg-paff=" << std::right << std::setw(5) << std::setprecision(2) << mean_placed_area_first_fail()
            << ", avg-pa=" << std::right << std::setw(5) << std::setprecision(2) << mean_placed_area()
            << ", med-pa=" << median_placed_area()
            << ", min-pa=" << min_placed_area()
            << ", max-pa=" << max_placed_area()
            << ", avg-pt=" << std::right << std::setw(5) << std::setprecision(2) << mean_placed_tiles()
            << ", med-pt=" << median_placed_tiles()
            << ", min-pt=" << min_placed_tiles()
            << ", max-pt=" << max_placed_tiles()
            << ", avg-fail=" << std::right << std::setw(8) << std::setprecision(2) << mean_failures()
            << ", max-fail=" << std::right << std::setw(5) << max_failures()
            << ", avg-alts=" << std::right << std::setw(5) << std::setprecision(2) << mean_alternatives()
            << ", max-alts=" << std::right << std::setw(3) << max_alternatives()
            << ")";
        return res.str();
    }
};

void report_csv(
        const std::map<std::pair<patchwork::PlacementPolicy, patchwork::PlacementEvaluation>, EvaluationStatistics> &stats);


int main(int argc, char **argv) {
    typedef std::pair<patchwork::PlacementPolicy, patchwork::PlacementEvaluation> ppp;
    std::map<ppp, EvaluationStatistics> stats;
    std::vector<ppp> keys;
    for (const auto &policy : patchwork::all_policies) {
        for (const auto &evaluation : patchwork::all_evaluations) {
            const auto key = ppp(policy, evaluation);
            keys.emplace_back(key);
            stats[key] = EvaluationStatistics(policy, evaluation);
        }
    }

    // Clock function used.
    auto now = [] { return std::chrono::steady_clock::now(); };

    const auto total_start = now();

    patchwork::Context global_context(42);
    auto *board = new patchwork::PatchworkBoard();
    board->status(global_context.status_statistics());

    std::vector<int> tiles;
    tiles.reserve(patchwork::tile_count());
    for (int tile = 0; tile < patchwork::tile_count(); ++tile) {
        tiles.emplace_back(tile);
    }

    auto rng = std::default_random_engine {}; // NOLINT(cert-msc32-c,cert-msc51-cpp)

    const int runs = 1000;

    for (int run = 0; run < runs; ++run) {
        for (const auto &key : keys) {
            const auto start = now();
            std::vector<int> alts;
            patchwork::Context context(run);
            auto *current = dynamic_cast<patchwork::PatchworkBoard *>(board->clone(context.clone_statistics()));
            int first_fail = tiles.size()+1;
            int placed_area_first_fail = 81;
            int tile_index = 0;
            for (const auto &tile : tiles) {
                auto placement_option = patchwork::place_tile(context, current, tile, key.first, key.second);
                ++tile_index;
                if (placement_option.has_value()) {
                    delete current;
                    current = placement_option.value().first;
                    alts.emplace_back(placement_option.value().second);
                    assert(!current->failed());
                    if (patchwork::DEBUG) {
                        current->print(std::cerr);
                    }
                } else {
                    first_fail = std::min(tile_index, first_fail);
                    placed_area_first_fail = std::min(placed_area_first_fail, current->placed_squares().min());
                }
            }
            const auto end = now();
            const std::chrono::duration<double, std::milli> duration =
                    end - start;
            int placed_squares = current->placed_squares().val();
            int placed_tiles = 0;
            for (const auto &tile : current->tiles()) {
                if (tile.assigned() && tile.val()) {
                    placed_tiles += 1;
                }
            }
            stats[key].add(duration.count(), first_fail, placed_area_first_fail,
                    placed_squares, placed_tiles, context.search_statistics().fail, alts);
            global_context += context;
        }
        std::shuffle(std::begin(tiles), std::end(tiles), rng);
        std::cerr << ".";
        if ((run != 0 && run % 100 == 0) || (run == runs-1)) {
            std::cerr << "  " << (run+1) << " of " << runs << " done" << std::endl;
        }
    }

    const auto total_end = now();
    const std::chrono::duration<double, std::milli> total_duration =
            total_end - total_start;
    std::cout.setf(std::ios::floatfield,std::ios::fixed);
    std::cout << "Running all " << stats.begin()->second.tests() << " tests on "
              << keys.size() << " combinations took " << total_duration.count() << " milliseconds." << std::endl;
    std::cout << "In total, " << global_context.search_statistics().node << " nodes were created, with "
              << global_context.search_statistics().fail << " failures and "
              << global_context.search_statistics().propagate << " propagations" << std::endl;

    for (const auto &stat : stats) {
        std::cout << stat.second.to_string() << std::endl;
    }

    report_csv(stats);

    return 0;
}

void report_csv(
        const std::map<std::pair<patchwork::PlacementPolicy, patchwork::PlacementEvaluation>, EvaluationStatistics> &stats) {
    using namespace std;
    using namespace patchwork;
    int precision = 10;

    cout << endl << endl
         << "************* Csv report BEGIN **************"
         << endl << endl;

    cout << "FullName,"
            "Placement,Evaluation,"
            "Tests,"
            "AvgRuntime,"
            "AvgRuntimePerAlternativ,"
            "AvgFirstFail,"
            "AvgFirstFailPlacedArea,"
            "AvgPlacedArea,"
            "MedianPlacedArea,"
            "MinPlacedArea,"
            "MaxPlacedArea,"
            "AvgPlacedTiles,"
            "MedianPlacedTiles,"
            "MinPlacedTiles,"
            "MaxPlacedTiles,"
            "AvgFail,"
            "MaxFail,"
            "AvgAlternatives,"
            "MaxAlternatives,"
         << endl;
    std::cout.setf(std::ios::floatfield,std::ios::fixed);
    for (const auto &stat : stats) {
        PlacementPolicy policy = stat.first.first;
        PlacementEvaluation evaluation = stat.first.second;
        EvaluationStatistics estat = stat.second;
        cout << to_string(policy) << "-" << to_string(evaluation) << ",";
        cout << to_string(policy) << "," << to_string(evaluation) << ",";
        cout << estat.tests() << ",";
        cout << setprecision(precision) << estat.mean_runtime() << ",";
        cout << setprecision(precision) << estat.mean_runtime_per_alt() << ",";
        cout << setprecision(precision) << estat.mean_first_fail() << ",";
        cout << setprecision(precision) << estat.mean_placed_area_first_fail() << ",";
        cout << setprecision(precision) << estat.mean_placed_area() << ",";
        cout << estat.median_placed_area() << ",";
        cout << estat.min_placed_area() << ",";
        cout << estat.max_placed_area() << ",";
        cout << setprecision(precision) << estat.mean_placed_tiles() << ",";
        cout << estat.median_placed_tiles() << ",";
        cout << estat.min_placed_tiles() << ",";
        cout << estat.max_placed_tiles() << ",";
        cout << setprecision(precision) << estat.mean_failures() << ",";
        cout << estat.max_failures() << ",";
        cout << setprecision(precision) << estat.mean_alternatives() << ",";
        cout << estat.max_alternatives() << ",";
        cout << endl;
    }

    cout << endl << endl
         << "************* Csv report END   **************"
         << endl << endl;
}
