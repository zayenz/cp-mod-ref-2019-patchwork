//
// Created by Mikael Zayenz Lagerkvist
//

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCUnusedGlobalDeclarationInspection"

#ifndef PATCHWORK_LIB_H
#define PATCHWORK_LIB_H

#include "tiles.h"
#include "board.h"

#include <gecode/driver.hh>
#include <gecode/int.hh>

#include <vector>
#include <cassert>
#include <optional>
#include <chrono>

namespace patchwork {
    static volatile bool DEBUG = false;

    class PatchworkBoard;

    /**
     * Modes define how to find the next state.
     *
     * The mode is embedded in the Context, enabling it to be threaded through the nested calls easily.
     */
    enum Mode {
        /// Find all possible next states
        MODE_ALL,
        /// Find a single next state, based on patch value and the current PlacementPolicy and PlacementEvaluation
        MODE_SINGLE,
        /// Count the amount of branching, using steps for income as the next state
        MODE_BRANCHING,
    };

    /**
     * The policy to use for placing a part.
     */
    enum PlacementPolicy {
        PP_NONE,
        PP_NONE_All_ROTATIONS,
        PP_FIRST_FAIL,
        PP_FIRST_FAIL_ALL_ROTATIONS,
        PP_AFC,
        PP_AFC_ALL_ROTATIONS,
        PP_ACTION,
        PP_ACTION_ALL_ROTATIONS,
        PP_CHB,
        PP_CHB_ALL_ROTATIONS,
        PP_AFC_SUM,
        PP_AFC_SUM_ALL_ROTATIONS,
        PP_ACTION_SUM,
        PP_ACTION_SUM_ALL_ROTATIONS,
        PP_CHB_SUM,
        PP_CHB_SUM_ALL_ROTATIONS,
        PP_BOTTOM_LEFT,
        PP_BOTTOM_LEFT_ALL_ROTATIONS,
        PP_BL_LB,
        PP_BL_LB_ALL_ROTATIONS,
        PP_PARETO_BOTTOM_LEFT,
        PP_PARETO_BOTTOM_LEFT_ALL_ROTATIONS,
        PP_ALL,
    };

    /**
     * The evaluation to use for determining which among many placements is the best placement.
     */
    enum PlacementEvaluation {
        PE_RND,
        PE_FIRST,
        PE_MIN_X_EXTENSION,
        PE_MIN_Y_EXTENSION,
        PE_MIN_EXTENSION,
        PE_MINIMUM_PLACEMENT_REGRET,
        PE_MAXIMUM_PLACEMENT_REGRET,
    };

    /**
     * A context contains a configuration and statistics-collection.
     */
    class Context {
    private:
        Gecode::CloneStatistics clone_statistics_;
        Gecode::CommitStatistics commit_statistics_;
        Gecode::Search::Statistics search_statistics_;
        Gecode::Search::Options search_options_;
        Gecode::Rnd rnd_;
        const Mode mode_;
        const PlacementPolicy placement_policy_;
        const PlacementEvaluation placement_evaluation_;
    public:
        explicit Context(
                int seed,
                Mode mode = MODE_ALL,
                PlacementPolicy policy = PP_ALL,
                PlacementEvaluation evaluation = PE_FIRST
        );

        Context with(Mode mode,
                     PlacementPolicy policy,
                     PlacementEvaluation evaluation
        );

        Gecode::CloneStatistics &clone_statistics();
        Gecode::StatusStatistics &status_statistics();
        Gecode::CommitStatistics &commit_statistics();
        Gecode::Search::Statistics &search_statistics();

        const Gecode::Search::Options &search_options() const;

        void reset_statistics() {
            clone_statistics_.reset();
            commit_statistics_.reset();
            search_statistics_.reset();
        }

        Context& operator+=(const Context& c) {
            clone_statistics_ += c.clone_statistics_;
            commit_statistics_ += c.commit_statistics_;
            search_statistics_ += c.search_statistics_;
            return *this;
        }

        Gecode::Rnd &rnd();

        const Mode mode() const;

        const PlacementPolicy policy() const;

        const PlacementEvaluation evaluation() const;

        void print_statistics(std::ostream &out, PatchworkBoard *root,
                              std::optional<std::chrono::duration<double, std::milli>> runtime = std::optional<std::chrono::duration<double, std::milli>>());
    };

    namespace {
        Gecode::StatusStatistics unused_status_statistics; // NOLINT(cert-err58-cpp)
        Gecode::CloneStatistics unused_clone_statistics;   // NOLINT(cert-err58-cpp)
    }

    /**
     * Gecode Space for representing a Patchwork 9 by 9 board to place patches/tiles on.
     */
    class PatchworkBoard : public Gecode::Space {
    private:
        /// Width of the board
        const int width_;
        /// Height of the board
        const int height_;
        /// Number of tiles to place
        const int ntiles_;
        /// Number of colors for board squares
        const int ncolors_;
        /// Number of board squares (including artificial end-of-line column)
        const int nsquares_;
        /// The value for empty squares
        const int empty_color_;
        /// The value for end squares
        const int end_color_;

        /// The variables for the board.
        Gecode::IntVarArray board_;

        /// The variables representing the chosen tiles.
        Gecode::BoolVarArray tile_is_used_;

        /// The variables representing the chosen tiles rotation.
        Gecode::IntVarArray tile_rotation_;

        /// Boolean variables representing placement of different tiles.
        /// The matrix is set up with columns representing individual squares in the board,
        /// and rows representing the colors/tiles.
        Gecode::BoolVarArray bool_board_;

        /// Boolean variables representing placement of any tile.
        /// True iff any element is placed. Smaller than board, no extra artificial column.
        Gecode::BoolVarArray monochrome_board_;

        /// The variables representing the usages of the tiles in columns and rows
        Gecode::IntVarArray usages_;

        /// The variables representing the first x-index/column of a tile (or one large if the part is not used)
        Gecode::IntVarArray first_x_;

        /// The variables representing the first y-index/row of a tile (or one large if the part is not used)
        Gecode::IntVarArray first_y_;

        /// The number of placed squares
        Gecode::IntVar placed_squares_;

        //
        // The following code for AFC/Action/CHB is inspired by the job-shop example in Gecode
        // https://github.com/Gecode/gecode/blob/master/examples/job-shop.cpp
        //

        /// AFC information
        Gecode::IntAFC board_afc_;
        /// AFC information for the bool board
        Gecode::BoolAFC bool_board_afc_;
        /// AFC-based cost
        double board_afc_over_size(int i) const {
            return board_[i].afc() / board_[i].size();
        }
        /// AFC-based cost, all variables
        double board_sum_afc_over_size(int i) const {
            double afc = board_[i].afc();
            for (int color = 0; color < ncolors_; ++color) {
                afc += bool_board_[i + color * ncolors_].afc();
            }
            return afc / board_[i].size();
        }
    public:
        /// Trampoline function for AFC-based cost
        static double afc_merit(const Space &home, const Gecode::BoolVar &, int i) {
            return dynamic_cast<const PatchworkBoard &>(home).board_afc_over_size(i);
        }
        /// Trampoline function for AFC-based cost, all variables
        static double afc_sum_merit(const Space &home, const Gecode::BoolVar &, int i) {
            return dynamic_cast<const PatchworkBoard &>(home).board_sum_afc_over_size(i);
        }
    private:
        /// Action information
        Gecode::IntAction board_action_;
        /// Action information for the bool board
        Gecode::BoolAction bool_board_action_;
        /// Action-based cost
        double board_action_over_size(int i) const {
            return board_action_[i] / board_[i].size();
        }
        /// Action-based cost, all variables
        double board_sum_action_over_size(int i) const {
            double action = board_action_[i];
            for (int color = 0; color < ncolors_; ++color) {
                action += bool_board_action_[i + color * ncolors_];
            }
            return action / board_[i].size();
        }
    public:
        /// Trampoline function for Action-based cost
        static double action_merit(const Space &home, const Gecode::BoolVar &, int i) {
            return dynamic_cast<const PatchworkBoard &>(home).board_action_over_size(i);
        }
        /// Trampoline function for Action-based cost, all variables
        static double action_sum_merit(const Space &home, const Gecode::BoolVar &, int i) {
            return dynamic_cast<const PatchworkBoard &>(home).board_sum_action_over_size(i);
        }
    private:
        /// CHB information
        Gecode::IntCHB board_chb_;
        /// CHB information
        Gecode::BoolCHB bool_board_chb_;
        /// CHB-based cost
        double board_chb_over_size(int i) const {
            return board_chb_[i] / board_[i].size();
        }
        /// CHB-based cost, all variables
        double board_sum_chb_over_size(int i) const {
            double chb = board_chb_[i];
            for (int color = 0; color < ncolors_; ++color) {
                chb += bool_board_chb_[i + color * ncolors_];
            }
            return chb / board_[i].size();
        }
    public:
        /// Trampoline function for CHB-based cost
        static double chb_merit(const Space &home, const Gecode::BoolVar &, int i) {
            return dynamic_cast<const PatchworkBoard &>(home).board_chb_over_size(i);
        }
        /// Trampoline function for CHB-based cost, all variables
        static double chb_sum_merit(const Space &home, const Gecode::BoolVar &, int i) {
            return dynamic_cast<const PatchworkBoard &>(home).board_sum_chb_over_size(i);
        }
    private:
        /// First fail
        double board_first_fail(int i) const {
            return board_[i].size();
        }
    public:
        /// Trampoline function for first fail
        static double first_fail_merit(const Space& home, const Gecode::BoolVar&, int i) {
            return dynamic_cast<const PatchworkBoard &>(home).board_first_fail(i);
        }
    public:
        /// Construction of the model.
        explicit PatchworkBoard();

        /// Constructor for cloning \a s
        PatchworkBoard(PatchworkBoard &s);

        /// Copy space during cloning
        PatchworkBoard *copy() override;

        /// Convenience method to create a clone with the correct return type
        PatchworkBoard *make_clone(Context& context);

        /// Convenience method to create a clone with the correct return type
        PatchworkBoard *make_clone(
                Gecode::StatusStatistics& status_statistics = unused_status_statistics,
                Gecode::CloneStatistics& clone_statistics = unused_clone_statistics
        );

        /// Print solution
        void print(std::ostream &os) const;

        int width() const { return width_; }
        int height() const { return height_; }
        int ntiles() const { return ntiles_; }
        int ncolors() const { return ncolors_; }
        int nsquares() const { return nsquares_; }
        int empty_color() const { return empty_color_; }
        int end_color() const { return end_color_; }
        const Gecode::IntVarArray &board() const { return board_; }
        const Gecode::BoolVarArray &tiles() const { return tile_is_used_; }
        const Gecode::IntVarArray &tile_rotation() const { return tile_rotation_; }
        const Gecode::BoolVarArray &bool_board() const { return bool_board_; }
        const Gecode::BoolVarArray &monochrome_board() const { return monochrome_board_; }
        const Gecode::IntVarArray &usages() const { return usages_; }
        const Gecode::IntVarArray &first_x() const { return first_x_; }
        const Gecode::IntVarArray &first_y() const { return first_y_; }
        const Gecode::IntVar &placed_squares() const { return placed_squares_; }


        int width() { return width_; }
        int height() { return height_; }
        int ntiles() { return ntiles_; }
        int ncolors() { return ncolors_; }
        int nsquares() { return nsquares_; }
        Gecode::IntVarArray &board() { return board_; }
        Gecode::BoolVarArray &tile_is_used() { return tile_is_used_; }
        Gecode::IntVarArray &tile_rotation() { return tile_rotation_; }
        Gecode::BoolVarArray &bool_board() { return bool_board_; }
        Gecode::BoolVarArray &monochrome_board() { return monochrome_board_; }
        Gecode::IntVarArray &usages() { return usages_; }
        Gecode::IntVarArray &first_x() { return first_x_; }
        Gecode::IntVarArray &first_y() { return first_y_; }
        Gecode::IntVar &placed_squares() { return placed_squares_; }

        const Gecode::IntVar& square(int x, int y) const {
            Gecode::Matrix<Gecode::IntVarArray> m(board_, width_, height_);
            return m(x, y);
        }

        Gecode::BoolVarArgs bool_board(int tile) {
            assert(valid_tile(tile));
            Gecode::Matrix<Gecode::BoolVarArray> m(bool_board_, nsquares_, ncolors_);
            Gecode::BoolVarArgs tile_board = m.row(tile);
            return tile_board;
        }

        bool valid_tile(int tile) const;

    };

    /**
     * Connect the variables so that first is the index fo the first non-zero variable in vars.
     *
     * If no variable in vars is non-zero, first should be the value of the index one past the last in vars.
     *
     * @param home the Home to post the constraint in
     * @param first The index of the first non-zero variable, or vars.size() is all in vars are zero
     * @param vars the variable array to inspect
     */
    void first_index(Gecode::Home home, const Gecode::IntVar& first, Gecode::IntVarArgs vars);

    /**
     * Prints a variable as a square
     *
     * @param os The stream to print to
     * @param square The variable representing the square to print
     */
    void print_square(std::ostream &os, const Gecode::IntVar &square);

    /**
     * @param policy The policy to name
     * @return The name of the policy
     */
    std::string to_string(PlacementPolicy policy);

    /**
     * @param evaluation The evaluation to name
     * @return The name of the evaluation
     */
    std::string to_string(PlacementEvaluation evaluation);

    const std::vector<PlacementPolicy> all_policies = { // NOLINT(cert-err58-cpp)
            PP_NONE,
            PP_NONE_All_ROTATIONS,
            PP_FIRST_FAIL,
            PP_FIRST_FAIL_ALL_ROTATIONS,
            PP_AFC,
            PP_AFC_ALL_ROTATIONS,
            PP_ACTION,
            PP_ACTION_ALL_ROTATIONS,
            PP_CHB,
            PP_CHB_ALL_ROTATIONS,
            PP_AFC_SUM,
            PP_AFC_SUM_ALL_ROTATIONS,
            PP_ACTION_SUM,
            PP_ACTION_SUM_ALL_ROTATIONS,
            PP_CHB_SUM,
            PP_CHB_SUM_ALL_ROTATIONS,
            PP_BOTTOM_LEFT,
            PP_BOTTOM_LEFT_ALL_ROTATIONS,
            PP_BL_LB,
            PP_BL_LB_ALL_ROTATIONS,
            PP_PARETO_BOTTOM_LEFT,
            PP_PARETO_BOTTOM_LEFT_ALL_ROTATIONS,
            PP_ALL,
    };

    const std::vector<PlacementEvaluation> all_evaluations = { // NOLINT(cert-err58-cpp)
            PE_RND,
            PE_FIRST,
            PE_MIN_X_EXTENSION,
            PE_MIN_Y_EXTENSION,
            PE_MIN_EXTENSION,
            PE_MINIMUM_PLACEMENT_REGRET,
            PE_MAXIMUM_PLACEMENT_REGRET,
    };

    /**
     *
     * @param start The start state to place the tile in (must not have any active branchings installed)
     * @param tile the tile to place
     * @param policy The policy to place the tile with
     * @return All placemenet in new clones derived from \a start with the tile \a tile placed according to the policy
     */
    std::optional<std::vector<PatchworkBoard *>>
    place_tile_all(Context &context, PatchworkBoard *start, int tile, PlacementPolicy policy);

    /**
     *
     * @param start The start state to place the tile in (must not have any active branchings installed)
     * @param tile the tile to place
     * @param policy The policy to place the tile with
     * @param evaluation The evaluation to use if the policy produces multiple choices
     * @return A new clone derived from \a start with the tile \a tile placed along with the number of alternatives
     *         choosen amon, or none if it is not possible
     */
    std::optional<std::pair<PatchworkBoard *, int>>
    place_tile(Context &context, PatchworkBoard *start, int tile, PlacementPolicy policy,
               PlacementEvaluation evaluation);

    enum Player {
        YELLOW = 0,
        GREEN = 1,
    };


    inline Player other(Player player) {
        switch (player) {
            case YELLOW:
                return GREEN;
            case GREEN:
                return YELLOW;
        }
    }

    static const int Players = 2;
    static const int TileCount = 33 + 5;
    static const int Squares = 9 * 9;

    class NextTiles;

    class State {
    private:
        friend NextTiles;

        bool finished_;
        PatchworkBoard* board_[Players];
        bool tiles_bought_[Players][TileCount];
        bool tiles_placed_[Players][TileCount];
        int bought_tile_count_[Players];
        int placed_tile_count_[Players];
        int bought_tile_area_[Players];
        int bought_tile_buttons_[Players];
        Player current_player_;
        int position_[Players];
        int purse_[Players];
        std::vector<int> tiles_;
        std::vector<int> singles_;
        int pos_;

        explicit State(Context& context) :
                finished_(false),
                board_{new PatchworkBoard(), new PatchworkBoard()},
                tiles_bought_{}, // Zero initialization of array
                tiles_placed_{}, // Zero initialization of array
                bought_tile_count_{0, 0},
                placed_tile_count_{0, 0},
                bought_tile_area_{0, 0},
                bought_tile_buttons_{0, 0},
                current_player_(YELLOW),
                position_{0, 0},
                purse_{5, 5},
                tiles_(patchwork::choosable_tiles()),
                singles_(patchwork::single_tiles()),
                pos_(1)
                {

            assert(TileCount == patchwork::tile_count());

            unsigned long size = tiles_.size();
            for (int i = 1; i < size-1; ++i) {
                const int j = i + context.rnd()(size - i);
                if (i != j) {
                    std::swap(tiles_[i], tiles_[j]);
                }
            }
            
            assert(patchwork::tile(tiles_[0]).is_start_tile());
        }

        State(Context& context, const State *state) :
                finished_(state->finished_),
                board_{state->board_[0]->make_clone(context),
                       state->board_[1]->make_clone(context)},
                tiles_bought_{}, // Zero initialization of array
                tiles_placed_{}, // Zero initialization of array
                bought_tile_count_{state->bought_tile_count_[0], state->bought_tile_count_[1]},
                placed_tile_count_{state->placed_tile_count_[0], state->placed_tile_count_[1]},
                bought_tile_area_{state->bought_tile_area_[0], state->bought_tile_area_[1]},
                bought_tile_buttons_{state->bought_tile_buttons_[0], state->bought_tile_buttons_[1]},
                current_player_(state->current_player_),
                position_{state->position_[0], state->position_[1]},
                purse_{state->purse_[0], state->purse_[1]},
                tiles_(state->tiles_),
                singles_(state->tiles_),
                pos_(state->pos_)
        {
            for (int player = 0; player < Players; ++player) {
                board_[player] = dynamic_cast<PatchworkBoard*>(state->board_[player]->clone());
                position_[player] = state->position_[player];
                purse_[player] = state->purse_[player];
                for (int tile = 0; tile < TileCount; ++tile) {
                    tiles_bought_[player][tile] = state->tiles_bought_[player][tile];
                    tiles_placed_[player][tile] = state->tiles_placed_[player][tile];
                }
            }
        }

        State(Context context, const State *state, const Player player, PatchworkBoard *board) :
                finished_(state->finished_),
                board_{nullptr, nullptr},
                tiles_bought_{}, // Zero initialization of array
                tiles_placed_{}, // Zero initialization of array
                bought_tile_count_{state->bought_tile_count_[0], state->bought_tile_count_[1]},
                placed_tile_count_{state->placed_tile_count_[0], state->placed_tile_count_[1]},
                bought_tile_area_{state->bought_tile_area_[0], state->bought_tile_area_[1]},
                bought_tile_buttons_{state->bought_tile_buttons_[0], state->bought_tile_buttons_[1]},
                current_player_(state->current_player_),
                position_{state->position_[0], state->position_[1]},
                purse_{state->purse_[0], state->purse_[1]},
                tiles_(state->tiles_),
                singles_(state->tiles_),
                pos_(state->pos_)
        {
            switch (player) {
                case YELLOW: {
                    board_[YELLOW] = board;
                    board_[GREEN] = state->board_[GREEN]->make_clone(context);
                    break;
                }
                case GREEN: {
                    board_[YELLOW] = state->board_[YELLOW]->make_clone(context);
                    board_[GREEN] = board;
                    break;
                }
            }

            for (int p = 0; p < Players; ++p) {
                board_[p] = dynamic_cast<PatchworkBoard*>(state->board_[p]->clone());
                position_[p] = state->position_[p];
                purse_[p] = state->purse_[p];
                for (int tile = 0; tile < TileCount; ++tile) {
                    tiles_bought_[p][tile] = state->tiles_bought_[p][tile];
                    tiles_placed_[p][tile] = state->tiles_placed_[p][tile];
                }
            }
        }



        void remove_tile(int pos, int tile) {
            assert(0 <= pos && pos < tiles_.size());
            assert(tiles_[pos] == tile);
            tiles_.erase(tiles_.begin() + pos);
        }

        std::optional<int> buy_tile(int pos, int tile) {
            std::optional<int> result;

            remove_tile(pos, tile);
            pos_ = pos % tiles_.size();
            const TileSource &source = patchwork::tile(tile);
            add_tile_to_current_player(tile, source);

            int old_pos = position_[current_player_];
            int new_pos = std::min(old_pos + source.time_cost(), patchwork::board_length-1);

            if (patchwork::board_collections[old_pos] != patchwork::board_collections[new_pos]) {
                // If the number of collectiosn changed between new and old pos, 
                // collect the buttons
                purse_[current_player_] += bought_tile_buttons_[current_player_];
            }

            if (patchwork::board_singles[old_pos] != patchwork::board_singles[new_pos] &&
                patchwork::board_singles[position_[other(current_player_)] != patchwork::board_singles[new_pos]]) {
                // if the number of singles left changes between the old and new pos, and the other player also has a 
                // different number of singles, get the next single and add to current players board
                assert(!singles_.empty());
                int single_tile = *singles_.begin();
                singles_.erase(singles_.begin());
                const TileSource &single_source = patchwork::tile(single_tile);
                add_tile_to_current_player(single_tile, single_source);
                result = std::make_optional(single_tile);
            }

            position_[current_player_] = new_pos;
            if (position_[current_player_] > position_[other(current_player_)]) {
                current_player_ = other(current_player_);
            }

            if (position_[0] == patchwork::board_length - 1 && position_[1] == patchwork::board_length - 1) {
                finished_ = true;
            }

            return result;
        }

        void add_tile_to_current_player(int tile, const TileSource &source) {
            tiles_bought_[current_player_][tile] = true;
            bought_tile_count_[current_player_] += 1;
            purse_[current_player_] -= source.button_cost();
            bought_tile_area_[current_player_] += source.area();
            bought_tile_buttons_[current_player_] += source.buttons();
            assert(purse_[current_player_] >= 0);
            PatchworkBoard *current_board = board_[current_player_];
            PatchworkBoard *other_board = board_[other(current_player_)];
            Gecode::rel(*current_board, current_board->tile_is_used()[tile], Gecode::IRT_EQ, 1);
            Gecode::rel(*other_board, other_board->tile_is_used()[tile], Gecode::IRT_EQ, 0);
        }

    public:
        virtual ~State();

        State *clone(Context& context) const {
            return new State(context, this);
        }

        static State *create(Context& context) {
            return new State(context);
        }

        NextTiles to_buy(Player player) const;

        bool is_placeable(Context& context, Player player, int tile) const;

        NextTiles select_placeable(Context& context, Player player, NextTiles next_tiles) const;

        Player current_player() const;

        std::vector<State *> next_states(Context &context) const {
            int alts = 0;
            return next_states(context, context.policy(), context.evaluation(), alts);
        }

        std::vector<State *> next_states(Context &context, int& alts) const {
            return next_states(context, context.policy(), context.evaluation(), alts);
        }

        std::vector<State *>
        next_states(Context &context,PlacementPolicy policy, PlacementEvaluation evaluation) const {
            int alts = 0;
            return next_states(context, policy, evaluation, alts);
        }

        std::vector<State *>
        next_states(Context &context,PlacementPolicy policy, PlacementEvaluation evaluation, int& alts) const;

        void step_for_buttons(Context &context);

        int steps_for_buttons_length() const;

        /**
         *
         * @param player
         * @param board The board to replace with. NOTE: TAKES OWNERSHIP OF THE ARGUMENT!!!!!!!
         * @return A new state, similar to the current one, with the player \a player's board replaced by \a board
         */
        State* with_board(Context& context, Player player, PatchworkBoard *board) const;

        /// Print solution
        void print(std::ostream &os) const;
        
        int current_score(Player player) const {
            return -2 * (Squares - bought_tile_area_[player]) + purse_[player];
        }

        int min_achievable_score(Player player) const {
            int remaining_income = bought_tile_buttons_[player] * (std::abs(
                    board_collections[board_length - 1] - board_collections[position_[player]]));
            return -2 * (Squares - bought_tile_area_[player]) + purse_[player] + remaining_income;
        }

        bool finished();
    };

    class NextTiles {
        std::optional<int> first_;
        std::optional<int> first_pos_;
        std::optional<int> second_;
        std::optional<int> second_pos_;
        std::optional<int> third_;
        std::optional<int> third_pos_;

        NextTiles() = default;
    public:
        explicit NextTiles(const State& state, int pos1 = -1, int pos2 = -1, int pos3 = -1) {
            if (pos1 >= 0) {
                first_ = std::optional<int>(state.tiles_[pos1]);
                first_pos_ = std::optional<int>(pos1);
            }
            if (pos2 >= 0) {
                assert(pos1 >= 0);
                second_ = std::optional<int>(state.tiles_[pos2]);
                second_pos_ = std::optional<int>(pos2);
            }
            if (pos3 >= 0) {
                assert(pos2 >= 0);
                third_ = std::optional<int>(state.tiles_[pos3]);
                third_pos_ = std::optional<int>(pos3);
            }
        }

        NextTiles(int tile1, int pos1, int tile2, int pos2, int tile3, int pos3) {
            if (pos1 >= 0) {
                first_ = std::optional<int>(tile1);
                first_pos_ = std::optional<int>(pos1);
            }
            if (pos2 >= 0) {
                assert(pos1 >= 0);
                second_ = std::optional<int>(tile2);
                second_pos_ = std::optional<int>(pos2);
            }
            if (pos3 >= 0) {
                assert(pos2 >= 0);
                third_ = std::optional<int>(tile3);
                third_pos_ = std::optional<int>(pos3);
            }
        }

        static NextTiles make_empty() {
            return NextTiles();
        }

        int size() const {
            if (third_.has_value()) {
                return 3;
            } else if (second_.has_value()) {
                return 2;
            } else if (first_.has_value()) {
                return 1;
            } else {
                return 0;
            }
        }

        bool empty() const {
            return !first_.has_value();
        }

        int tile1() const {
            return first_.value();
        }

        int tile2() const {
            return second_.value();
        }

        int tile3() const {
            return third_.value();
        }

        int pos1() const {
            return first_pos_.value();
        }

        int pos2() const {
            return second_pos_.value();
        }

        int pos3() const {
            return third_pos_.value();
        }

        int pos(int i) const {
            switch (i) {
                case 0:
                    return pos1();
                case 1:
                    return pos2();
                case 2:
                    return pos3();
                default:
                    assert(false);
                    GECODE_NEVER
                    return -1;
            }
        }

        int tile(int i) const {
            switch (i) {
                case 0:
                    return tile1();
                case 1:
                    return tile2();
                case 2:
                    return tile3();
                default:
                    assert(false);
                    GECODE_NEVER
                    return -1;
            }
        }
    };

}

#endif //PATCHWORK_LIB_H

#pragma clang diagnostic pop