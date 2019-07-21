//
// Created by Mikael Zayenz Lagerkvist
//

#include "lib.h"
#include "tiles.h"
#include "board.h"

#include <gecode/driver.hh>
#include <gecode/int.hh>
#include <gecode/minimodel.hh>

#include <iostream>
#include <iomanip>
#include <utility>

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCUnusedGlobalDeclarationInspection"
#pragma ide diagnostic ignored "OCDFAInspection"

using namespace Gecode;

namespace patchwork {
    PatchworkBoard::PatchworkBoard()
            : Gecode::Space(),
              width_(9 + 1), // Add one for extra row at end.
              height_(9),
              ntiles_(patchwork::tile_count()),
              ncolors_(ntiles_ + 2),
              nsquares_(width_ * height_),
              empty_color_(ntiles_),
              end_color_(ntiles_ + 1),
              board_(*this, nsquares_, 0, ncolors_),
              tile_is_used_(*this, ntiles_, 0, 1),
              tile_rotation_(*this, ntiles_, 0, 8),
              bool_board_(*this, ncolors_ * nsquares_, 0, 1),
              monochrome_board_(*this, (width_ - 1) * height_, 0, 1),
              usages_(*this, ntiles_ * (width_ + height_), 0, std::max(width_ - 1, height_)),
              first_x_(*this, ntiles_, 0, 9 + 1),
              first_y_(*this, ntiles_, 0, 9 + 1),
              placed_squares_(*this, 0, (width_ - 1) * height_),
              board_afc_(*this, board_, 0.99),
              bool_board_afc_(*this, bool_board_, 0.99),
              board_action_(*this, board_, 0.99),
              bool_board_action_(*this, bool_board_, 0.99),
              board_chb_(*this, board_),
              bool_board_chb_(*this, bool_board_)
    {

        // Matrix for the board (includes the end of line column)
        Matrix<IntVarArray> b(board_, width_, height_);

        // Matrix for the monochrome board (does not include the end of line column)
        Matrix<BoolVarArray> mb(monochrome_board_, width_ - 1, height_);

        // Boolean variables representing placement of different tiles.
        // The matrix is set up with columns representing individual squares in the board,
        // and rows representing the board for individual colors/tiles.
        Matrix<BoolVarArray> bool_board_matrix(bool_board_, nsquares_, ncolors_);

        // Matrix for the column and row usages for each tile.
        // Column c of the matrix represents the tile c
        Matrix<IntVarArray> usage_matrix(usages_, ntiles_, width_ + height_);


        // Set end-of-line markers
        for (int h = 0; h < height_; ++h) {
            for (int w = 0; w < width_ - 1; ++w) {
                rel(*this, b(w, h), IRT_NQ, end_color_);
            }
            rel(*this, b(width_ - 1, h), IRT_EQ, end_color_);
        }


        // Post channeling constraints
        //

        // Connect monochrome_board and board
        assert(end_color_ >= ntiles_ && empty_color_ >= ntiles_);
        for (int x = 0; x < width_ - 1; ++x) {
            for (int y = 0; y < height_; ++y) {
                rel(*this, b(x, y), IRT_LE, ntiles_, mb(x, y));
            }
        }

        // Connect the board and the bool board.
        // Each column represents a square, with one element in the column being true at the index of the corresponding board square variables value.
        for (int i = 0; i < nsquares_; ++i) {
            // All Boolean variables for a square
            BoolVarArgs square = bool_board_matrix.col(i);
            channel(*this, square, board_[i]);
        }

        // Connect the usages with the bool board
        for (int tile = 0; tile < ntiles_; ++tile) {
            BoolVarArgs tile_board = bool_board_matrix.row(tile);
            Matrix<BoolVarArgs> tile_board_matrix(tile_board, width_, height_);
            IntVarArgs tile_usages = usage_matrix.col(tile);
            int pos = 0;
            for (int col = 0; col < width_; ++col) {
                linear(*this, tile_board_matrix.col(col), IRT_EQ, tile_usages[pos]);
                ++pos;
            }
            for (int row = 0; row < height_; ++row) {
                linear(*this, tile_board_matrix.row(row), IRT_EQ, tile_usages[pos]);
                ++pos;
            }
        }

        // Connect first x and y with usages
        for (int tile = 0; tile < ntiles_; ++tile) {
            IntVarArgs tile_usages = usage_matrix.col(tile); // 9 column usages, 1 empty marker, and 9 row usages
            IntVarArgs col_usages = tile_usages.slice(0, 1, 9); // The first 9 elements, the column usages
            IntVarArgs row_usages = tile_usages.slice(10); // The last 9 elements, the row usages
            assert(row_usages.size() == 9);
            first_index(*this, first_x_[tile], col_usages);
            first_index(*this, first_y_[tile], row_usages);
        }

        // Post placement related constraints
        //

        // Set up reified rotation variables for use in placement and usage constraints
        // Each of the ntiles_ rows of the matrix corresponds to the 9 Booleans for a rotation
        // 0 for the tile not being used, 1-8 for the possible rotations
        BoolVarArgs all_rotation(*this, ntiles_ * 9, 0, 1);
        Matrix<BoolVarArgs> all_rotation_matrix(all_rotation, 9, ntiles_);
        for (int tile = 0; tile < ntiles_; ++tile) {
            // One row per tile, with each row consisting of 9 Booleans for the possible rotations
            BoolVarArgs rotation = all_rotation_matrix.row(tile);
            // Connect with the int var for the rotation
            channel(*this, rotation, tile_rotation_[tile]);
            // Connect with the control variables for if the tile is used or not
            // Alternative 0 is true iff the tile is not used.
            rel(*this, rotation[0], IRT_NQ, tile_is_used_[tile]);
        }

        // Post reified placement expressions for each tile
        for (int tile = 0; tile < ntiles_; ++tile) {
            // Boolean variables for the board representing if a part is placed
            BoolVarArgs tile_board = bool_board_matrix.row(tile);
            const REG &placement = patchwork::tile(tile).as_placement_expression();
            // Make placement expression reified with an additional first control variable and no-alternative rotation
            REG reified_placement =
                    (REG(1) + // Placement for control variable
                     placement // Placement on board
                    ) |
                    (REG(0) + // No placement for control variable
                     REG(1) + REG(0)(8, 8) + // "Alternative" 0/no alternative
                     REG(0)(nsquares_, nsquares_) // No placement on board
                    );
            const BoolVarArgs reified_tile_variables = tile_is_used_[tile] + all_rotation_matrix.row(tile) + tile_board;
            extensional(*this, reified_tile_variables, reified_placement);
        }

        // Post reified usage constraint for each tile
        for (int tile = 0; tile < ntiles_; ++tile) {
            const REG &usage = patchwork::tile(tile).usage_expression();
            // Make usage expression reified with an additional first control variable and no-alternative rotation
            REG reified_usage =
                    (REG(1) + // Placement for control variable
                     usage // Placement on board
                    ) |
                    (REG(0) + // No placement for control variable
                     REG(1) + REG(0)(8, 8) + // "Alternative" 0/no alternative
                     REG(0)(width_ + height_, width_ + height_) // No placement on board
                    );
            const IntVarArgs tile_usage = usage_matrix.col(tile);
            const BoolVarArgs bool_rotation = all_rotation_matrix.row(tile);
            IntVarArgs bool_rotation_as_int;
            for (const auto &var : bool_rotation) {
                bool_rotation_as_int << channel(*this, var);
            }
            const IntVarArgs reified_usage_variables =
                    channel(*this, tile_is_used_[tile]) + bool_rotation_as_int + tile_usage;
            channel(*this, tile_is_used_[tile]) + bool_rotation_as_int + tile_usage;
            extensional(*this, reified_usage_variables, reified_usage);
        }


        // Post tile placement sum both from reified variables and from the squares
        IntArgs tile_sizes;
        for (int tile = 0; tile < ntiles_; ++tile) {
            tile_sizes << patchwork::tile(tile).area();
        }
        linear(*this, tile_sizes, tile_is_used_, IRT_EQ, placed_squares_);
        linear(*this, monochrome_board_, IRT_EQ, placed_squares_);
    }

    PatchworkBoard::PatchworkBoard(PatchworkBoard &s) :
            Space(s), width_(s.width_), height_(s.height_),
            ntiles_(s.ntiles_), ncolors_(s.ncolors_), nsquares_(s.nsquares_),
            empty_color_(s.empty_color_), end_color_(s.end_color_),
            board_afc_(s.board_afc_), bool_board_afc_(s.bool_board_afc_),
            board_action_(s.board_action_), bool_board_action_(s.bool_board_action_),
            board_chb_(s.board_chb_), bool_board_chb_(s.bool_board_chb_)
    {
        board_.update(*this, s.board_);
        tile_is_used_.update(*this, s.tile_is_used_);
        tile_rotation_.update(*this, s.tile_rotation_);
        bool_board_.update(*this, s.bool_board_);
        monochrome_board_.update(*this, s.monochrome_board_);
        usages_.update(*this, s.usages_);
        first_x_.update(*this, s.first_x_);
        first_y_.update(*this, s.first_y_);
        placed_squares_.update(*this, s.placed_squares_);
    }

    PatchworkBoard *PatchworkBoard::copy() {
        return new PatchworkBoard(*this);
    }

    PatchworkBoard *PatchworkBoard::make_clone(Context &context) {
        return make_clone(context.status_statistics(), context.clone_statistics());
    }

    PatchworkBoard *PatchworkBoard::make_clone(
            Gecode::StatusStatistics &status_statistics,
            Gecode::CloneStatistics &clone_statistics) {
        if (status(status_statistics) == SS_FAILED) {
            return nullptr;
        }
        return dynamic_cast<PatchworkBoard *>(clone(clone_statistics));
    }

    void PatchworkBoard::print(std::ostream &os) const {
        for (int h = 0; h < height_; ++h) {
            os << "\t";
            for (int w = 0; w < width_ - 1; ++w) {
                print_square(os, board_[h * width_ + w]);
            }
            os << std::endl;
        }
        os << std::endl;
        os << "Tiles used : " << tile_is_used_ << std::endl;
        os << "First x : " << first_x_ << std::endl;
    }

    bool PatchworkBoard::valid_tile(int tile) const { return 0 <= tile && tile < ntiles_; }

    void print_square(std::ostream &os, const IntVar &square) {
        if (square.assigned()) {
            int val = square.val();
            if (val < patchwork::tile_count()) {
                char c = static_cast<char>(val < 10 ? '0' + val : 'A' + (val - 10));
                os << c;
            } else {
                os << ' ';
            }
        } else {
            os << "·";
        }
    }


    void first_index(Home home, const IntVar& first, IntVarArgs vars) {
        // The following decomposition is somewhat similar to the one in the global constraints catalog
        // http://sofdem.github.io/gccat/gccat/Clength_first_sequence.html
        assert(vars.size() > 0);

        // Set up variables indicating if the value in vars is a zero or not,
        BoolVarArgs is_zero(home, vars.size(), 0, 1);
        for (int i = 0; i < vars.size(); ++i) {
            rel(home, vars[i], IRT_EQ, 0, is_zero[i]);
        }

        // Set up variables indicating if the index is a part of the first run of zeroes
        BoolVarArgs is_first_run(home, vars.size(), 0, 1);
        rel(home, is_first_run[0], IRT_EQ, is_zero[0]); // Position 0 in run iff the first element is a zero
        for (int i = 1; i < vars.size(); ++i) {
            // An index is part of the first run iff the current element is a zero and the previous
            // index was part of the first run.
            rel(home, is_first_run[i] == (is_zero[i] && is_first_run[i - 1]));
        }

        // The first k variables in is_first_run will be true (and the rest false) when the first non-zero variable
        // in vars is at index k. When there is no such index, (that is, vars i zero everywhere) is_first_run
        // is true everywhere corresponding to the desired sentinel value for no index.
        // Thus a simple sum suffices to connect is_first_run and first
        linear(home, is_first_run, IRT_EQ, first);
    }


    /**
     * Find the first solution according to dfs search order.
     *
     * @param context The context for the operations (out parameter)
     * @param root The root space to start search from
     * @return A pointer to the solution, if any
     */
    PatchworkBoard *dfs(Context &context, PatchworkBoard *root);

    /**
     * Find the first solution according to dfs search order.
     *
     * @param context The context for the operations (out parameter)
     * @param root The root space to start search from
     * @param result Vector in which to store all the results (out parameter)
     */
    void dfs(Context &context, PatchworkBoard *root,
            /* out parameter, result collection */ std::vector<PatchworkBoard *> &result);

    /**
     * Find all solutions according to dfs search order.
     *
     * @param context The context for the operations (out parameter)
     * @param root The root space to start search from
     * @param result Vector in which to store all the results (out parameter)
     */
    void dfs_all(Context &context, PatchworkBoard *root,
            /* out parameter, result collection */ std::vector<PatchworkBoard *> &result);

    /**
     * For all possible rotations for tile \a tile, find the first solution according to dfs
     * search order.
     *
     * @param context The context for the operations (out parameter)
     * @param root The root space to start search from
     * @param tile The tile to do the search for
     * @param result Vector in which to store all the results (out parameter)
     */
    void dfs_all_rotations(Context &context, PatchworkBoard *root, int tile,
            /* out parameter, result collection */ std::vector<PatchworkBoard *> &result);

    /**
     *
     * @param root The space to investigave
     * @return The maximum value for x for which a variable is assigned.
     */
    int max_assigned_x_value(PatchworkBoard *root);

    /**
     *
     * @param root The space to investigave
     * @return The maximum value for y for which a variable is assigned.
     */
    int max_assigned_y_value(PatchworkBoard *root);

    /**
     * Using a policy, place the tile on the board.
     *
     * @param context  The context for the operations (out parameter)
     * @param root The root space to start placements in
     * @param tile The tile to place
     * @param policy The policy to use for placements
     * @param tile_board The board variables in \a root for the tile
     * @return A list of placement results
     */
    std::vector<PatchworkBoard *>
    place_tile_with_policy(Context &context, PatchworkBoard *root, int tile, PlacementPolicy policy,
                             const BoolVarArgs &tile_board);

    /**
     *
     * @param root The original space the placement were made in
     * @param tile The placed tile
     * @param placements The placements made
     * @param evaluation THe evaluation to use
     * @return The list of placements, re-ordered according to the \a evaluation
     */
    std::vector<PatchworkBoard *>
    evaluate_placement(Context context, PatchworkBoard *root, int tile, std::vector<PatchworkBoard *> placements,
                       PlacementEvaluation evaluation);

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCSimplifyInspection"

    std::optional<std::vector<PatchworkBoard *>>
    place_tile_all(Context &context, PatchworkBoard *start, int tile,
                   PlacementPolicy policy) {
        assert(start->valid_tile(tile) && "pre: tile must be a valid tile index");

        std::optional<std::vector<PatchworkBoard *>> no_solution;

        switch (start->status(context.status_statistics())) {
            case SS_FAILED:
                // Since the space is failed, it is not possible to place the tile
                return no_solution;
            case SS_SOLVED:
                // This is the desired outcome
                break;
            case SS_BRANCH:
                assert(false && "Pre: A root space should not include branchers");
                // The return here is just to ensure a reasonable result when asserts are turned off.
                return no_solution;
        }

        auto *root = start->make_clone(context);

        BoolVar &tile_used = root->tile_is_used()[tile];
        const BoolVarArgs &tile_board = root->bool_board(tile);
        if (tile_used.assigned()) {
            if (tile_used.val() == false) {
                // The tile is not placeable
                delete root;
                return no_solution;
            }
        } else {
            rel(*root, tile_used, IRT_EQ, true);
            if (root->status(context.status_statistics()) == SS_FAILED) {
                context.search_statistics().fail += 1;
                delete root;
                return no_solution;
            }
        }

        if (tile_board.assigned()) {
            // The tile is used and all the squares for the tile board are assigned.
            // Placement is already done!
            return std::make_optional(std::vector<PatchworkBoard*>{root});
        }


        std::vector<PatchworkBoard *> tile_placements = place_tile_with_policy(context, root, tile, policy,
                                                                                 tile_board);
        bool root_in_results = false;
        for (PatchworkBoard* tile_placement : tile_placements) {
            if (tile_placement == root) {
                root_in_results = true;
                break;
            }
        }
        if (!root_in_results) {
            delete root;
        }

        return std::make_optional(tile_placements);
    }
    
    std::optional<std::pair<PatchworkBoard *, int>>
    place_tile(Context &context, PatchworkBoard *start, int tile,
               PlacementPolicy policy, PlacementEvaluation evaluation) {
        assert(start->valid_tile(tile) && "pre: tile must be a valid tile index");

        std::optional<std::pair<PatchworkBoard *, int>> no_solution;

        switch (start->status(context.status_statistics())) {
            case SS_FAILED:
                // Since the space is failed, it is not possible to place the tile
                return no_solution;
            case SS_SOLVED:
                // This is the desired outcome
                break;
            case SS_BRANCH:
                assert(false && "Pre: A root space should not include branchers");
                // The return here is just to ensure a reasonable result when asserts are turned off.
                return no_solution;
        }

        auto *root = start->make_clone(context);

        BoolVar &tile_used = root->tile_is_used()[tile];
        const BoolVarArgs &tile_board = root->bool_board(tile);
        if (tile_used.assigned()) {
            if (tile_used.val() == false) {
                // The tile is not placeable
                delete root;
                return no_solution;
            }
        } else {
            rel(*root, tile_used, IRT_EQ, true);
            if (root->status(context.status_statistics()) == SS_FAILED) {
                context.search_statistics().fail += 1;
                delete root;
                return no_solution;
            }
        }

        if (tile_board.assigned()) {
            // The tile is used and all the squares for the tile board are assigned.
            // Placement is already done!
            return std::make_optional(std::make_pair(root, 1));
        }


        std::vector<PatchworkBoard *> tile_placements = place_tile_with_policy(context, root, tile, policy,
                                                                                 tile_board);
        if (tile_placements.empty()) {
            return no_solution;
        }

        std::vector<PatchworkBoard *> ordered_tile_placements = evaluate_placement(context, root, tile, tile_placements,
                                                                                   evaluation);
        assert(tile_placements.size() == ordered_tile_placements.size());

        bool root_in_placements = false;
        for (PatchworkBoard *placement : tile_placements) {
            if (placement == root) {
                root_in_placements = true;
                break;
            }
        }
        if (!root_in_placements) {
            delete root;
        }

        const std::optional<std::pair<PatchworkBoard *, int>> result =
                std::make_optional(std::make_pair(ordered_tile_placements[0], ordered_tile_placements.size()));

        for (int i = 1; i < ordered_tile_placements.size(); ++i) {
            delete ordered_tile_placements[i];
        }

        return result;
    }

    std::string to_string(PlacementPolicy policy) {
        switch (policy) {
            case PP_NONE:
                return "None";
            case PP_NONE_All_ROTATIONS:
                return "None(rot)";
            case PP_FIRST_FAIL:
                return "FirstFail";
            case PP_FIRST_FAIL_ALL_ROTATIONS:
                return "FirstFail(rot)";
            case PP_AFC:
                return "AFC/Size";
            case PP_AFC_ALL_ROTATIONS:
                return "AFC/size(rot)";
            case PP_ACTION:
                return "Action/Size";
            case PP_ACTION_ALL_ROTATIONS:
                return "Action/Size(rot)";
            case PP_CHB:
                return "CHB/Size";
            case PP_CHB_ALL_ROTATIONS:
                return "CHB/size(rot)";
            case PP_AFC_SUM:
                return "AFCSum/Size";
            case PP_AFC_SUM_ALL_ROTATIONS:
                return "AFCSum/size(rot)";
            case PP_ACTION_SUM:
                return "ActionSum/Size";
            case PP_ACTION_SUM_ALL_ROTATIONS:
                return "ActionSum/Size(rot)";
            case PP_CHB_SUM:
                return "CHBSum/Size";
            case PP_CHB_SUM_ALL_ROTATIONS:
                return "CHBSum/size(rot)";
            case PP_BOTTOM_LEFT:
                return "BottomLeft";
            case PP_BOTTOM_LEFT_ALL_ROTATIONS:
                return "BottomLeft(rot)";
            case PP_BL_LB:
                return "BL/LB";
            case PP_BL_LB_ALL_ROTATIONS:
                return "BL/LB(rot)";
            case PP_PARETO_BOTTOM_LEFT:
                return "ParetoBottomLeft";
            case PP_PARETO_BOTTOM_LEFT_ALL_ROTATIONS:
                return "ParetoBottomLeft(rot)";
            case PP_ALL:
                return "All";
        }
    }

    std::string to_string(PlacementEvaluation evaluation) {
        switch (evaluation) {
            case PE_RND:
                return "Random";
            case PE_FIRST:
                return "First";
            case PE_MIN_X_EXTENSION:
                return "MinXExtension";
            case PE_MIN_Y_EXTENSION:
                return "MinYExtension";
            case PE_MIN_EXTENSION:
                return "MinXYAreaExtension";
            case PE_MINIMUM_PLACEMENT_REGRET:
                return "PlacementRegret";
            case PE_MAXIMUM_PLACEMENT_REGRET:
                return "ReversePlacementRegret";
        }
    }

    std::vector<PatchworkBoard *>
    evaluate_placement(Context context, PatchworkBoard *root, int tile, std::vector<PatchworkBoard *> placements,
                       PlacementEvaluation evaluation) {
        if (placements.size() <= 1) {
            return placements;
        }

        using namespace std;
        int root_x = max_assigned_x_value(root);
        int root_y = max_assigned_y_value(root);
        typedef pair<PatchworkBoard *, int> ppi;
        vector<ppi> pwc;
        pwc.reserve(placements.size());
        for (int i = 0; i < placements.size(); ++i) {
            int value;
            switch (evaluation) {
                case PE_RND:
                    // This corresponds to no evaluation, just use input order
                    value = context.rnd()(2 * placements.size());
                    break;
                case PE_FIRST:
                    // This corresponds to no evaluation, just use input order
                    value = i;
                    break;
                case PE_MIN_X_EXTENSION: {
                    // the minimum increase of the assigned variables x extension
                    int x = max_assigned_x_value(placements[i]);
                    value = x - root_x;
                    break;
                }
                case PE_MIN_Y_EXTENSION: {
                    // the minimum increase of the assigned variables x extension
                    int y = max_assigned_y_value(placements[i]);
                    value = y - root_y;
                    break;
                }
                case PE_MIN_EXTENSION: {
                    // the minimum increase of the area of the bounding box of the assigned variables
                    int x = max_assigned_x_value(placements[i]);
                    int y = max_assigned_y_value(placements[i]);
                    value = (x * y) - root_x * root_y;
                    break;
                }
                case PE_MINIMUM_PLACEMENT_REGRET: {
                    // The minimum reduction in possibilities for board placements
                    int regret = 0;
                    for (int square_index = 0; square_index < root->board().size(); ++square_index) {
                        IntVar &root_square = root->board()[square_index];
                        IntVar &square = placements[i]->board()[square_index];

                        bool assigned_in_root = root_square.assigned();
                        bool assigned = square.assigned();
                        bool is_tile_placement = assigned && square.val() == tile;

                        if (!assigned_in_root && !is_tile_placement) {
                            int local_regret = static_cast<long long>(root_square.size()) -
                                               static_cast<long long>(square.size());
                            assert(local_regret >= 0 && "inv: propagation should never add values...");
                            regret += local_regret;
                        }
                    }
                    value = regret;
                    break;
                }
                case PE_MAXIMUM_PLACEMENT_REGRET: {
                    // The minimum reduction in possibilities for board placements
                    int regret = 0;
                    for (int square_index = 0; square_index < root->board().size(); ++square_index) {
                        IntVar &root_square = root->board()[square_index];
                        IntVar &square = placements[i]->board()[square_index];

                        bool assigned_in_root = root_square.assigned();
                        bool assigned = square.assigned();
                        bool is_tile_placement = assigned && square.val() == tile;

                        if (!assigned_in_root && !is_tile_placement) {
                            int local_regret = static_cast<long long>(root_square.size()) -
                                               static_cast<long long>(square.size());
                            assert(local_regret >= 0 && "inv: propagation should never add values...");
                            regret += local_regret;
                        }
                    }
                    value = -regret;
                    break;
                }
            }
            pwc.emplace_back(make_pair(placements[i], value));
        }
        stable_sort(pwc.begin(), pwc.end(), [](const ppi &a, const ppi &b) {
            return a.second < b.second;
        });

        vector<PatchworkBoard *> result;
        result.reserve(placements.size());
        for (const auto &item : pwc) {
            result.emplace_back(item.first);
        }

        if (DEBUG) {
            for (int i = 0; i < pwc.size(); ++i) {
                std::cerr << "i:" << std::right << std::setw(2) << i
                          << " s:" << std::right << std::setw(3) << pwc[i].second
                          << "   ";
            }
            std::cerr << std::endl;
            for (int y = 0; y < root->height(); ++y) {
                for (int i = 0; i < result.size(); ++i) {
                    for (int x = 0; x < root->width(); ++x) {
                        print_square(std::cerr, result[i]->square(x, y));
                    }
                    if (i < result.size() - 1) {
                        std::cerr << "   ";
                    } else {
                        std::cerr << std::endl;
                    }
                }
            }
            std::cerr << std::endl;
        }

        return result;
    }

    std::vector<PatchworkBoard *>
    place_tile_with_policy(Context &context, PatchworkBoard *root, int tile, PlacementPolicy policy,
                             const BoolVarArgs &tile_board) {
        BrancherGroup brancher_group;
        const Home bghome = Home(*root)(brancher_group);

        std::vector<PatchworkBoard *> result;
        switch (policy) {
            case PP_NONE: {
                // Add branching that just assigns variables at the start to contain the tile
                branch(bghome, tile_board, BOOL_VAR_NONE(), BOOL_VAL_MAX());
                dfs(context, root, result);
                break;
            }
            case PP_NONE_All_ROTATIONS: {
                // Add branching that just assigns variables at the start to contain the tile
                branch(bghome, tile_board, BOOL_VAR_NONE(), BOOL_VAL_MAX());
                dfs_all_rotations(context, root, tile, result);
                break;
            }
            case PP_FIRST_FAIL: {
                // Add branching that uses standard first-fail on the board variable size to assign
                // the tile_board variables
                branch(bghome, tile_board, BOOL_VAR_MERIT_MAX(&PatchworkBoard::first_fail_merit), BOOL_VAL_MAX());
                dfs(context, root, result);
                break;
            }
            case PP_FIRST_FAIL_ALL_ROTATIONS: {
                // Add branching that uses standard first-fail on the board variable size to assign
                // the tile_board variables
                branch(bghome, tile_board, BOOL_VAR_MERIT_MAX(&PatchworkBoard::first_fail_merit), BOOL_VAL_MAX());
                dfs_all_rotations(context, root, tile, result);
                break;
            }
            case PP_AFC: {
                // Add branching that uses standard AFC over domain size on the board variable size to assign
                // the tile_board variables
                branch(bghome, tile_board, BOOL_VAR_MERIT_MAX(&PatchworkBoard::afc_merit), BOOL_VAL_MAX());
                dfs(context, root, result);
                break;
            }
            case PP_AFC_ALL_ROTATIONS: {
                // Add branching that uses standard AFC over domain size on the board variable size to assign
                // the tile_board variables
                branch(bghome, tile_board, BOOL_VAR_MERIT_MAX(&PatchworkBoard::afc_merit), BOOL_VAL_MAX());
                dfs_all_rotations(context, root, tile, result);
                break;
            }
            case PP_ACTION: {
                // Add branching that uses standard action over domain size  on the board variable size to assign
                // the tile_board variables
                branch(bghome, tile_board, BOOL_VAR_MERIT_MAX(&PatchworkBoard::action_merit), BOOL_VAL_MAX());
                dfs(context, root, result);
                break;
            }
            case PP_ACTION_ALL_ROTATIONS: {
                // Add branching that uses standard action over domain size  on the board variable size to assign
                // the tile_board variables
                branch(bghome, tile_board, BOOL_VAR_MERIT_MAX(&PatchworkBoard::action_merit), BOOL_VAL_MAX());
                dfs_all_rotations(context, root, tile, result);
                break;
            }
            case PP_CHB: {
                // Add branching that uses standard chb over domain size  on the board variable size to assign
                // the tile_board variables
                branch(bghome, tile_board, BOOL_VAR_MERIT_MAX(&PatchworkBoard::chb_merit), BOOL_VAL_MAX());
                dfs(context, root, result);
                break;
            }
            case PP_CHB_ALL_ROTATIONS: {
                // Add branching that uses standard chb  over domain size  on the board variable size to assign
                // the tile_board variables
                branch(bghome, tile_board, BOOL_VAR_MERIT_MAX(&PatchworkBoard::chb_merit), BOOL_VAL_MAX());
                dfs_all_rotations(context, root, tile, result);
                break;
            }
            case PP_AFC_SUM: {
                // Add branching that uses standard AFC over domain size on the board variable size to assign
                // the tile_board variables
                branch(bghome, tile_board, BOOL_VAR_MERIT_MAX(&PatchworkBoard::afc_sum_merit), BOOL_VAL_MAX());
                dfs(context, root, result);
                break;
            }
            case PP_AFC_SUM_ALL_ROTATIONS: {
                // Add branching that uses standard AFC over domain size on the board variable size to assign
                // the tile_board variables
                branch(bghome, tile_board, BOOL_VAR_MERIT_MAX(&PatchworkBoard::afc_sum_merit), BOOL_VAL_MAX());
                dfs_all_rotations(context, root, tile, result);
                break;
            }
            case PP_ACTION_SUM: {
                // Add branching that uses standard action over domain size  on the board variable size to assign
                // the tile_board variables
                branch(bghome, tile_board, BOOL_VAR_MERIT_MAX(&PatchworkBoard::action_sum_merit), BOOL_VAL_MAX());
                dfs(context, root, result);
                break;
            }
            case PP_ACTION_SUM_ALL_ROTATIONS: {
                // Add branching that uses standard action over domain size  on the board variable size to assign
                // the tile_board variables
                branch(bghome, tile_board, BOOL_VAR_MERIT_MAX(&PatchworkBoard::action_sum_merit), BOOL_VAL_MAX());
                dfs_all_rotations(context, root, tile, result);
                break;
            }
            case PP_CHB_SUM: {
                // Add branching that uses standard chb over domain size  on the board variable size to assign
                // the tile_board variables
                branch(bghome, tile_board, BOOL_VAR_MERIT_MAX(&PatchworkBoard::chb_sum_merit), BOOL_VAL_MAX());
                dfs(context, root, result);
                break;
            }
            case PP_CHB_SUM_ALL_ROTATIONS: {
                // Add branching that uses standard chb  over domain size  on the board variable size to assign
                // the tile_board variables
                branch(bghome, tile_board, BOOL_VAR_MERIT_MAX(&PatchworkBoard::chb_sum_merit), BOOL_VAL_MAX());
                dfs_all_rotations(context, root, tile, result);
                break;
            }
            case PP_BOTTOM_LEFT: {
                // The bottom left heuristic here uses the fact that a minimizing branching order
                // will when combined with dfs exploration find the "optimal"/minimal assignment as the first solution

                // Add branching first on x position, then the y position, and finally on the values
                branch(bghome, root->first_x()[tile], INT_VAL_MIN());
                branch(bghome, root->first_y()[tile], INT_VAL_MIN());
                // Finish up by assigning all the tile board values
                branch(bghome, tile_board, BOOL_VAR_NONE(), BOOL_VAL_MAX());
                dfs(context, root, result);
                break;
            }
            case PP_BOTTOM_LEFT_ALL_ROTATIONS: {
                // The bottom left heuristic here uses the fact that a minimizing branching order
                // will when combined with dfs exploration find the "optimal"/minimal assignment as the first solution

                // Add branching first on x position, then the y position, and finally on the values
                branch(bghome, root->first_x()[tile], INT_VAL_MIN());
                branch(bghome, root->first_y()[tile], INT_VAL_MIN());
                // Finish up by assigning all the tile board values
                branch(bghome, tile_board, BOOL_VAR_NONE(), BOOL_VAL_MAX());
                dfs_all_rotations(context, root, tile, result);
                break;
            }
            case PP_BL_LB: {
                // The bottom left/left bottom  heuristic here uses the fact that a minimizing branching order
                // will when combined with dfs exploration find the "optimal"/minimal assignment as the first solution
                {
                    // Add branching first on x position, then the y position, and finally on the values
                    branch(bghome, root->first_x()[tile], INT_VAL_MIN());
                    branch(bghome, root->first_y()[tile], INT_VAL_MIN());
                    // Finish up by assigning all the tile board values
                    branch(bghome, tile_board, BOOL_VAR_NONE(), BOOL_VAL_MAX());
                    dfs(context, root, result);
                    brancher_group.kill(*root);
                }
                {
                    // Add branching first on y position, then the x position, and finally on the values
                    branch(bghome, root->first_y()[tile], INT_VAL_MIN());
                    branch(bghome, root->first_x()[tile], INT_VAL_MIN());
                    // Finish up by assigning all the tile board values
                    branch(bghome, tile_board, BOOL_VAR_NONE(), BOOL_VAL_MAX());
                    dfs(context, root, result);
                    brancher_group.kill(*root);
                }
                break;
            }
            case PP_BL_LB_ALL_ROTATIONS: {
                // The bottom left/left bottom  heuristic here uses the fact that a minimizing branching order
                // will when combined with dfs exploration find the "optimal"/minimal assignment as the first solution
                {
                    // Add branching first on x position, then the y position, and finally on the values
                    branch(bghome, root->first_x()[tile], INT_VAL_MIN());
                    branch(bghome, root->first_y()[tile], INT_VAL_MIN());
                    // Finish up by assigning all the tile board values
                    branch(bghome, tile_board, BOOL_VAR_NONE(), BOOL_VAL_MAX());
                    dfs_all_rotations(context, root, tile, result);
                    brancher_group.kill(*root);
                }
                {
                    // Add branching first on y position, then the x position, and finally on the values
                    branch(bghome, root->first_y()[tile], INT_VAL_MIN());
                    branch(bghome, root->first_x()[tile], INT_VAL_MIN());
                    // Finish up by assigning all the tile board values
                    branch(bghome, tile_board, BOOL_VAR_NONE(), BOOL_VAL_MAX());
                    dfs_all_rotations(context, root, tile, result);
                    brancher_group.kill(*root);
                }
                break;
            }
            case PP_PARETO_BOTTOM_LEFT: {
                // The bottom left heuristic here uses the fact that a minimizing branching order will when combined with dfs exploration
                // find the "optimal"/minimal assignment as the first solution

                // This creates the pareto front of all the bottom left positions, by iterating over all possible x-values manually,
                // and finding the minimum y-position for all rotations
                int max_assigned_x = max_assigned_x_value(root);
                // Add branching first the y position. THe x position will be assigned manually
                branch(bghome, root->first_y()[tile], INT_VAL_MIN());
                // Finish up by assigning all the tile board values
                branch(bghome, tile_board, BOOL_VAR_NONE(), BOOL_VAL_MAX());
                for (int x = 0; x <= std::min(max_assigned_x + 1, root->width() - 1); ++x) {
                    auto *clone = dynamic_cast<PatchworkBoard *>(root->clone(context.clone_statistics()));
                    if (clone == nullptr) {
                        assert(false && "inv: cloning should work");
                        continue;
                    }
                    rel(*clone, clone->first_x()[tile], IRT_EQ, x);
                    if (clone->status(context.status_statistics()) == SS_FAILED) {
                        context.search_statistics().fail += 1;
                        // Assigning to the x value did not work.
                        delete clone;
                        continue;
                    }
                    dfs(context, clone, result);
                    delete clone;
                }
                break;
            }
            case PP_PARETO_BOTTOM_LEFT_ALL_ROTATIONS: {
                // The bottom left heuristic here uses the fact that a minimizing branching order will when combined with dfs exploration
                // find the "optimal"/minimal assignment as the first solution

                // This creates the pareto front of all the bottom left positions, by iterating over all possible x-values manually,
                // and finding the minimum y-position for all rotations
                int max_assigned_x = max_assigned_x_value(root);
                // Add branching first the y position. THe x position will be assigned manually
                branch(bghome, root->first_y()[tile], INT_VAL_MIN());
                // Finish up by assigning all the tile board values
                branch(bghome, tile_board, BOOL_VAR_NONE(), BOOL_VAL_MAX());
                for (int x = 0; x <= std::min(max_assigned_x + 1, root->width() - 1); ++x) {
                    auto *clone = dynamic_cast<PatchworkBoard *>(root->clone(context.clone_statistics()));
                    if (clone == nullptr) {
                        assert(false && "inv: cloning should work");
                        continue;
                    }
                    rel(*clone, clone->first_x()[tile], IRT_EQ, x);
                    if (clone->status(context.status_statistics()) == SS_FAILED) {
                        context.search_statistics().fail += 1;
                        // Assigning to the x value did not work.
                        delete clone;
                        continue;
                    }
                    dfs_all_rotations(context, clone, tile, result);
                    delete clone;
                }
                break;
            }
            case PP_ALL: {
                // Add branching that just assigns variables at the start to contain the tile
                branch(bghome, tile_board, BOOL_VAR_NONE(), BOOL_VAL_MAX());
                dfs_all(context, root, result);
                break;
            }
        }

        // Clean up all installed branchers
        brancher_group.kill(*root);

        return result;
    }


    int max_assigned_x_value(PatchworkBoard *root) {
        int max_assigned_x = 0;
        Matrix<IntVarArray> board(root->board(), root->width(), root->height());
        for (int x = 0; x < root->width(); ++x) {
            for (int y = 0; y < root->height(); ++y) {
                if (board(x, y).assigned()) {
                    max_assigned_x = x;
                    break;
                }
            }
            if (max_assigned_x < x) {
                break;
            }
        }
        return max_assigned_x;
    }

    int max_assigned_y_value(PatchworkBoard *root) {
        int max_assigned_y = 0;
        Matrix<IntVarArray> board(root->board(), root->width(), root->height());
        for (int y = 0; y < root->height(); ++y) {
            for (int x = 0; x < root->width(); ++x) {
                if (board(x, y).assigned()) {
                    max_assigned_y = y;
                    break;
                }
            }
            if (max_assigned_y < y) {
                break;
            }
        }
        return max_assigned_y;
    }

    void dfs(Context &context, PatchworkBoard *root, std::vector<PatchworkBoard *> &result) {
        DFS<PatchworkBoard> engine(root, context.search_options());
        PatchworkBoard *placement = engine.next();
        if (placement != nullptr) {
            if (placement->status(context.status_statistics()) != SS_FAILED) {
                result.emplace_back(placement);
            } else {
                context.search_statistics().fail += 1;
                delete placement;
            }
        }
        context.search_statistics() += engine.statistics();
    }

    PatchworkBoard *dfs(Context &context, PatchworkBoard *root) {
        DFS<PatchworkBoard> engine(root, context.search_options());
        PatchworkBoard *placement = engine.next();
        context.search_statistics() += engine.statistics();
        if (placement != nullptr) {
            if (placement->status(context.status_statistics()) != SS_FAILED) {
                return placement;
            } else {
                context.search_statistics().fail += 1;
                delete placement;
                return nullptr;
            }
        } else {
            return nullptr;
        }
    }

    void dfs_all(Context &context, PatchworkBoard *root, std::vector<PatchworkBoard *> &result) {
        DFS<PatchworkBoard> engine(root, context.search_options());
        PatchworkBoard *placement = engine.next();
        while (placement != nullptr) {
            if (placement->status(context.status_statistics()) != SS_FAILED) {
                result.emplace_back(placement);
            } else {
                context.search_statistics().fail += 1;
                delete placement;
            }
            placement = engine.next();
        }
        context.search_statistics() += engine.statistics();
    }

    void dfs_all_rotations(Context &context, PatchworkBoard *root, int tile, std::vector<PatchworkBoard *> &result) {
        for (IntVarValues i(root->tile_rotation()[tile]); i(); ++i) {
            auto *clone = dynamic_cast<PatchworkBoard *>(root->clone(context.clone_statistics()));
            if (clone == nullptr) {
                assert(false && "inv: cloning should work");
                continue;
            }
            rel(*clone, clone->tile_rotation()[tile], IRT_EQ, i.val());
            if (clone->status(context.status_statistics()) == SS_FAILED) {
                context.search_statistics().fail += 1;
                // Assigning to the rotation value did not work.
                delete clone;
                continue;
            }
            DFS<PatchworkBoard> engine(clone, context.search_options());
            PatchworkBoard *placement = engine.next();
            delete clone;
            if (placement != nullptr) {
                if (placement->status(context.status_statistics()) != SS_FAILED) {
                    result.emplace_back(placement);
                } else {
                    context.search_statistics().fail += 1;
                    delete placement;
                }
            }
            context.search_statistics() += engine.statistics();
        }
    }

#pragma clang diagnostic pop

    CloneStatistics &Context::clone_statistics() {
        return clone_statistics_;
    }

    StatusStatistics &Context::status_statistics() {
        return search_statistics_;
    }

    CommitStatistics &Context::commit_statistics() {
        return commit_statistics_;
    }

    const Search::Options &Context::search_options() const {
        return search_options_;
    }

    Context::Context(
            const int seed,
            Mode mode,
            PlacementPolicy policy,
            PlacementEvaluation evaluation
    )
            : rnd_(seed), mode_(mode),
            placement_policy_(policy), placement_evaluation_(evaluation) {
        search_options_.threads = 1;
    }

    Search::Statistics &Context::search_statistics() {
        return search_statistics_;
    }

    void Context::print_statistics(std::ostream &out, PatchworkBoard *root,
                                   std::optional<std::chrono::duration<double, std::milli>> runtime) {
        unsigned int n_p = PropagatorGroup::all.size(*root);
        unsigned int n_b = BrancherGroup::all.size(*root);

        out << "Initial" << std::endl
            << "\tpropagators: " << n_p << std::endl
            << "\tbranchers:   " << n_b << std::endl
            << std::endl
            << "Summary" << std::endl;
        if (runtime.has_value()) {
            out << "\truntime (ms): " << runtime.value().count() << std::endl;
        }
        out << "\tpropagations: " << search_statistics_.propagate << std::endl
            << "\tnodes:        " << search_statistics_.node << std::endl
            << "\tfailures:     " << search_statistics_.fail << std::endl
            << "\trestarts:     " << search_statistics_.restart << std::endl
            << "\tno-goods:     " << search_statistics_.nogood << std::endl
            << "\tpeak depth:   " << search_statistics_.depth << std::endl;
    }

    Rnd &Context::rnd() {
        return rnd_;
    }
    const Mode Context::mode() const {
        return mode_;
    }

    const PlacementPolicy Context::policy() const {
        return placement_policy_;
    }

    const PlacementEvaluation Context::evaluation() const {
        return placement_evaluation_;
    }

    Context Context::with(Mode mode, PlacementPolicy policy, PlacementEvaluation evaluation) {
        return Context(rnd_.seed(), mode, policy, evaluation);
    }

    NextTiles State::to_buy(Player player) const {
        if (tiles_.empty()) {
            return NextTiles(*this);
        }

        int pos1 = pos_;
        int pos2 = tiles_.size() > 1 ? ((pos_ + 1) % tiles_.size()) : -1;
        int pos3 = tiles_.size() > 2 ? ((pos_ + 2) % tiles_.size()) : -1;

        if (pos3 != -1) {
            const TileSource &source3 = patchwork::tile(tiles_[pos3]);
            if (source3.button_cost() > purse_[player]) {
                pos3 = -1;
            }
        }

        if (pos2 != -1) {
            const TileSource &source2 = patchwork::tile(tiles_[pos2]);
            if (source2.button_cost() > purse_[player]) {
                pos2 = pos3;
                pos3 = -1;
            }
        }

        const TileSource &source1 = patchwork::tile(tiles_[pos1]);
        if (source1.button_cost() > purse_[player]) {
            pos1 = pos2;
            pos2 = pos3;
            pos3 = -1;
        }

        return NextTiles(*this, pos1, pos2, pos3);
    }

    bool State::is_placeable(Context &context, Player player, int tile) const {
        if (tiles_placed_[player][tile]) {
            // Already placed
            return true;
        }
        if (tiles_bought_[other(player)][tile]) {
            // Implies tile placed implies tile bought
            // Other player has the tile, not placeable
            return false;
        }
        const BoolVar &tile_used = board_[player]->tile_is_used()[tile];
        if (tile_used.assigned() && !tile_used.val()) {
            // Propagation says the tile can't be used
            return false;
        }

        if (bought_tile_count_[player] == placed_tile_count_[player]) {
            // Fast mode, just check result of propagation
            // All bought tiles are placed, so propagation will set unplaceable tiles to false
            return tile_used.max() == 1;
        } else {
            // Slow mode, search for answer
            // Clone the state, set the tile to used and all undecided tiles to not used
            // Find solution for the board variables
            PatchworkBoard *clone = board_[player]->make_clone(context);
            for (int t = 0; t < clone->ntiles(); ++t) {
                BoolVar &clone_tile_is_used = clone->tile_is_used()[t];
                if (t == tile) {
                    rel(*clone, clone_tile_is_used, IRT_EQ, 1);
                } else if (!clone_tile_is_used.assigned()) {
                    rel(*clone, clone_tile_is_used, IRT_EQ, 0);
                }
            }
            assert(clone->tile_is_used().assigned() && "All tiles must have assigned state for usage here!");
            // Tests on placement model indicate that min domain size (first fail) is one of (if not the) fastest
            // ways to make placements.
            branch(*clone, clone->board(), INT_VAR_SIZE_MIN(), INT_VAL_MAX());

            if (clone->status(context.status_statistics()) == SS_FAILED) {
                delete clone;
                return false;
            }

            PatchworkBoard *placement = dfs(context, clone);

            const bool found_placement = placement != nullptr;

            delete clone;
            delete placement; // nullptr deletion is ok

            return found_placement;
        }
    }

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCSimplifyInspection"

    NextTiles State::select_placeable(Context &context, Player player, NextTiles next_tiles) const {
        int pos1 = -1;
        int pos2 = -1;
        int pos3 = -1;
        int tile1 = -1;
        int tile2 = -1;
        int tile3 = -1;
        switch (next_tiles.size()) {
            case 3:
                if (is_placeable(context, player, next_tiles.tile3())) {
                    pos3 = next_tiles.pos3();
                    tile3 = next_tiles.tile3();
                }
                // Fall through intended
            case 2:
                if (is_placeable(context, player, next_tiles.tile2())) {
                    pos2 = next_tiles.pos2();
                    tile2 = next_tiles.tile2();
                }
                // Fall through intended
            case 1:
                if (is_placeable(context, player, next_tiles.tile1())) {
                    pos1 = next_tiles.pos1();
                    tile1 = next_tiles.tile1();
                }
                break;
            case 0:
                break;
            default:
                assert(false && "only 0 to 3 valid sizes");
        }

        if (pos2 == -1) {
            pos2 = pos3;
            tile2 = tile3;
            pos3 = -1;
            tile3 = -1;
        }
        if (pos1 == -1) {
            pos1 = pos2;
            tile1 = tile2;
            pos2 = pos3;
            tile2 = tile3;
            pos3 = -1;
            tile3 = -1;
        }

        return NextTiles(tile1, pos1, tile2, pos2, tile3, pos3);
    }

    std::vector<State *>
    State::next_states(Context &context, PlacementPolicy policy, PlacementEvaluation evaluation, int& alts) const {
        std::vector<State *> next_states;

        const Player current_player = current_player_; // Remember the player we are making moves from
        // In particular, buy_tile will update who the current player is

        // Create all the moves.
        const NextTiles &tiles = select_placeable(context, current_player, to_buy(current_player));

        switch (context.mode()) {
            case MODE_ALL: {
                next_states.reserve(1 + 60 * tiles.size());
                // Create step forward move
                State *stepped = clone(context);
                stepped->step_for_buttons(context);
                next_states.emplace_back(stepped);

                for (int i = 0; i < tiles.size(); ++i) {
                    State *next = clone(context);
                    auto tile = tiles.tile(i);
                    (void) next->buy_tile(tiles.pos(i), tile);
                    auto all_placements = place_tile_all(context, next->board_[current_player], tile, policy);
                    if (all_placements.has_value()) {
                        for (auto placement : all_placements.value()) {
                            next_states.emplace_back(next->with_board(context, current_player, placement));
                        }
                    }
                }

                break;
            }

            case MODE_SINGLE: {
                next_states.reserve(1);
                alts = 1;
                const int step_income = steps_for_buttons_length();
                int choice = -1; // We let negative represent steps forward as choice 
                double best_value = step_income;
                for (int i = 0; i < tiles.size(); ++i) {
                    double value = patchwork::tile(tiles.tile(i)).value(position_[current_player]);
                    if (value > best_value) {
                        best_value = value;
                        choice = i;
                    }
                }
                State *next = clone(context);
                if (choice == -1) {
                    next->step_for_buttons(context);
                    next_states.emplace_back(next);
                } else {
                    int tile = tiles.tile(choice);
                    const std::optional<int> single_tile_gained = next->buy_tile(tiles.pos(choice), tile);

                    const auto &placed = place_tile(context, next->board_[current_player], tile,
                                                    policy, evaluation);
                    delete next->board_[current_player];
                    assert(placed.has_value());
                    auto[board, _alts] = placed.value();
                    alts += _alts;
                    next->tiles_placed_[current_player][tile] = true;
                    next->placed_tile_count_[current_player] += 1;
                    if (single_tile_gained.has_value()) {
                        int single_tile = single_tile_gained.value();
                        auto placed_single = place_tile(context, board, single_tile, context.policy(),
                                                        context.evaluation());
                        delete board;
                        auto[board_with_single, _alts_for_single] = placed_single.value();
                        alts += _alts_for_single;
                        next->tiles_placed_[current_player][single_tile] = true;
                        next->placed_tile_count_[current_player] += 1;
                        next->board_[current_player] = board_with_single;
                    } else {
                        next->board_[current_player] = board;
                    }

                    next_states.emplace_back(next);
                }
                break;
            }

            case MODE_BRANCHING: {
                alts = 1;

                for (int choice = 0; choice < tiles.size(); ++choice) {
                    State *test = clone(context);

                    int tile = tiles.tile(choice);
                    const std::optional<int> single_tile_gained = test->buy_tile(tiles.pos(choice), tile);
                    const auto &placed = place_tile(context, test->board_[current_player], tile,
                                                    policy, evaluation);
                    auto[board, _alts] = placed.value();
                    delete test->board_[current_player];
                    assert(placed.has_value());
                    alts = _alts;
                    test->tiles_placed_[current_player][tile] = true;
                    test->placed_tile_count_[current_player] += 1;
                    if (single_tile_gained.has_value()) {
                        int single_tile = single_tile_gained.value();
                        auto placed_single = place_tile(context, board, single_tile, context.policy(),
                                                        context.evaluation());
                        delete board;
                        auto[board_with_single, _alts_for_single] = placed_single.value();
                        alts += _alts_for_single;
                        test->tiles_placed_[current_player][single_tile] = true;
                        test->placed_tile_count_[current_player] += 1;
                        test->board_[current_player] = board_with_single;
                    } else {
                        test->board_[current_player] = board;
                    }

                    delete test;
                }

                State *result = clone(context);
                result->step_for_buttons(context);
                next_states.emplace_back(result);

                break;
            }

        }

        return next_states;
    }

    void State::step_for_buttons(Context &context) {
        int forward = steps_for_buttons_length();
        purse_[current_player_] += forward;
        position_[current_player_] += forward;
        current_player_ = other(current_player_);
        if (position_[0] == board_length - 1 && position_[1] == board_length - 1) {
            finished_ = true;
        }
    }

    int State::steps_for_buttons_length() const {
        int cpos = position_[current_player_];
        int opos = position_[other(current_player_)];
        assert(cpos <= opos);
        int forward = opos - cpos + 1;
        if (opos == board_length-1) {
            forward -= 1;
        }
        assert(forward > 0);
        return forward;
    }

    State::~State() {
        delete board_[0];
        delete board_[1];

    }

    State* State::with_board(Context& context, const Player player, PatchworkBoard *board) const {
        return new State(context, this, player, board);
    }

    void State::print(std::ostream &os) const {
        auto ps = [](Player player) {
            return player == GREEN ? "G" : "Y";
        };
        auto fill = [&](Player player, int step) {
            return position_[player] == step ? ps(player) : " ";
        };
        os << (finished_ ? "Game finished." : "Game in progress.")
           << " Current player is " << ps(current_player_) << std::endl;
        for (Player p : {GREEN, YELLOW}) {
            os << ps(p) << ": [" << fill(p, 0) << "]";
            for (int step = 1; step < board_length; ++step) {
                if (board_singles[step - 1] != board_singles[step]) {
                    if (position_[GREEN] < step && position_[YELLOW] < step) {
                        os << "#";
                    } else {
                        os << " ";
                    }
                }
                if (board_collections[step - 1] != board_collections[step]) {
                        os << "B";
                }
                os << "[" << fill(p, step) << "]";
            }
            os << std::endl;
        }

        for (Player p : {GREEN, YELLOW}) {
            os << "Stats for " << ps(p) << ": ";
            os << "Score: " << current_score(p) << " / " << min_achievable_score(p);
            os << ". Tile area: " << bought_tile_area_[p] << ", count: " << bought_tile_count_[p] << ", buttons: " << bought_tile_buttons_[p];
            os << ". Tiles: ";
            for (int tile = 0; tile < tile_count(); ++tile) {
                if (tiles_bought_[p][tile]) {
                    os << tile << " ";
                }
            }
            os << std::endl;
        }
        os << " Board Green  Board Yellow" << std::endl;
        for (int row = 0; row < 9; ++row) {
            os << "  ";
            for (Player p : {GREEN, YELLOW}) {
                for (int col = 0; col < 9; ++col) {
                    print_square(os, board_[p]->square(col, row));
                }
                os << "    ";
            }
            os << std::endl;
        }

    }

    bool State::finished() {
        return finished_;
    }

    Player State::current_player() const {
        return current_player_;
    }
}

#pragma clang diagnostic pop