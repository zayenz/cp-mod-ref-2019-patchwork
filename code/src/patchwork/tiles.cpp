//
// Created by Mikael Zayenz Lagerkvist
//

#include <iostream>
#include <gecode/minimodel.hh>
#include <gecode/int.hh>
#include <gecode/driver.hh>
#include <utility>

#include "lib.h"
#include "tiles.h"
#include "symmetry.h"
#include "board.h"

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCUnusedGlobalDeclarationInspection"

namespace patchwork {
    Tile::Tile(int width, int height, std::vector<bool> marks) : width(width), height(height), marks(std::move(marks)) {
        assert(this->marks.size() == width * height);
    }

    const Gecode::REG Tile::make_placement_expression(const Gecode::REG& mark, const Gecode::REG& empty_mark,
            const int board_width, const int board_height) const {
        using Gecode::REG;
        REG empty = empty_mark;

        // The fixed separation between rows, including the extra end-of-line column
        const auto fixed_separation_length = board_width + 1 - width;
        const REG fixed_separation = empty(fixed_separation_length, fixed_separation_length);

        // Start anywhere (that is, with arbitrary number of empty squares)
        REG result = *empty;

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                if (at(x, y)) {
                    result += mark;
                } else {
                    result += empty;
                }
            }
            if (y < height - 1) {
                // Between rows, add the fixed separation
                result += fixed_separation;
            }
        }

        result += *empty;

        return result;
    }

    const Gecode::REG Tile::make_usage_expression() const {
        using Gecode::REG;
        REG empty(0);

        REG result;

        // Add column usages
        result += *empty;
        for (const auto &col_usage : compute_column_usage()) {
            result += REG(col_usage);
        }
        result += *empty;

        // Add separation marker
        result += empty;

        // Add row usages
        result += *empty;
        for (const auto &row_usage : compute_row_usage()) {
            result += REG(row_usage);
        }
        result += *empty;

        return result;
    }

    const std::vector<int> Tile::compute_column_usage() const noexcept {
        std::vector<int> usage(width);

        for (int col = 0; col < width; ++col) {
            int sum = 0;
            for (int row = 0; row < height; ++row) {
                sum += at(col, row);
            }
            usage[col] = sum;
        }

        return usage;
    }

    const std::vector<int> Tile::compute_row_usage() const noexcept {
        std::vector<int> usage(height);

        for (int row = 0; row < height; ++row) {
            int sum = 0;
                for (int col = 0; col < width; ++col) {
                sum += at(col, row);
            }
            usage[row] = sum;
        }

        return usage;
    }

    TileSource::TileSource(const int id,
                           const int width, const int height, const int amount,
                           const int button_cost, const int time_cost, const int buttons,
                           const char *tile_pattern)
            : id_(id),
              amount_(amount),
              area_(count_area(width, height, tile_pattern)),
              button_cost_(button_cost),
              time_cost_(time_cost),
              buttons_(buttons),
              alternatives_(make_unique_tiles(width, height, tile_pattern)),
              placement_expression_(make_placement_expression(alternatives_)),
              usage_expression_(make_usage_expression(alternatives_))
              {}

    std::vector<Tile> TileSource::make_unique_tiles(const int width, const int height, const char *tile_pattern) {
        assert(width * height == strlen(tile_pattern));

        const std::vector<symmetry::vsymmfunc> symmetries {
            symmetry::id, symmetry::rot90, symmetry::rot180, symmetry::rot270,
            symmetry::flipx, symmetry::flipy, symmetry::flipd1, symmetry::flipd2
        };

        std::vector<bool> base_tile;
        base_tile.reserve(width * height);
        for (int i = 0; i < width * height; ++i) {
            base_tile.emplace_back(tile_pattern[i] == 'X');
        }

        std::vector<Tile> result;
        result.reserve(symmetries.size());

        for (const auto &symmetry : symmetries) {
            std::vector<bool> transformed(width * height);
            int transformed_width, transformed_height;
            symmetry(base_tile, width, height,
                    transformed, transformed_width, transformed_height);

            Tile tile(transformed_width, transformed_height, transformed);

            bool unique = true;
            for (const auto &alternative : result) {
                if (tile == alternative) {
                    unique = false;
                    break;
                }
            }
            if (unique) {
                result.emplace_back(tile);
            }
        }

        return result;
    }

    const int TileSource::count_area(const int width, const int height, const char *pattern) {
        int result = 0;
        for (int i = 0; i < width * height; ++i) {
            if (pattern[i] == 'X') {
                ++result;
            }
        }
        return result;
    }



    const Gecode::REG TileSource::make_placement_expression(const std::vector<Tile>& alternatives) {
        using Gecode::REG;
        const REG mark(1);
        const REG empty(0);

        REG result;

        for (int i = 0; i < alternatives.size(); ++i) {
            REG alternative_index = make_alternative_index(i);
            REG placement_expression = alternatives[i].make_placement_expression(mark, empty, 9, 9);
            result = result | (alternative_index + placement_expression);
        }
        
        return result;
    }

    const Gecode::REG TileSource::make_usage_expression(const std::vector<Tile>& alternatives) {
        using Gecode::REG;
        const REG mark(1);
        const REG empty(0);

        REG result;

        for (int i = 0; i < alternatives.size(); ++i) {
            REG alternative_index = make_alternative_index(i);
            REG placement_expression = alternatives[i].make_usage_expression();
            result = result | (alternative_index + placement_expression);
        }

        return result;
    }

    const Gecode::REG TileSource::make_alternative_index(int alternative) {
        assert(0 <= alternative && alternative < 8);
        using Gecode::REG;
        REG yes(1);
        REG no(0);

        REG result;
        result += no; // Use 1-based indexing
        for (int i = 0; i < 8; ++i) {
            if (i == alternative) {
                result += yes;
            } else {
                result += no;
            }
        }

        return result;
    }
    
    const Gecode::REG TileSource::as_placement_expression() const {
        return placement_expression_;
    }

    const Gecode::REG TileSource::usage_expression() const {
        return usage_expression_;
    }

    const int TileSource::id() const {
        return id_;
    }

    const int TileSource::amount() const {
        return amount_;
    }

    const int TileSource::area() const {
        return area_;
    }

    const int TileSource::button_cost() const {
        return button_cost_;
    }

    const int TileSource::time_cost() const {
        return time_cost_;
    }

    const int TileSource::buttons() const {
        return buttons_;
    }

    const bool TileSource::is_start_tile() const {
        return area() == 2;
    }

    const bool TileSource::is_bonus_tile() const {
        return area() == 1;
    }

    const double TileSource::value(int square) const {
        assert(0 <= square && square < board_length - 1);
        const int remaining_time = board_length - square;
        const int remaining_collections = board_collections[board_length-1] - board_collections[square];
        double score = static_cast<double>(2 * area() - button_cost() + remaining_collections)
                / std::min(time_cost(), remaining_time);
        return score;
    }


    const std::vector<TileSource> base_tiles = { // NOLINT(cert-err58-cpp)
            {0,
             1, 1, 5,
             0, 0, 0,
             "X"
            },
            {1,
             2, 1, 1,
             2, 1, 0,
             "XX"
            },
            {2,
             2, 2, 1,
             3, 1, 0,
             " X"
             "XX"
            },
            {3,
             2, 2, 1,
             1, 3, 0,
             " X"
             "XX"
            },
            {4,
             2, 2, 1,
             6, 5, 2,
             " X"
             "XX"
            },
            {5,
             3, 1, 1,
             2, 2, 0,
             "XXX"
            },
            {6,
             3, 2, 1,
             3, 2, 1,
             " XX"
             "XX "
            },
            {7,
             3, 2, 1,
             7, 6, 3,
             " XX"
             "XX "
            },
            {8,
             3, 2, 1,
             4, 6, 2,
             "  X"
             "XXX"
            },
            {9,
             3, 2, 1,
             4, 2, 1,
             "  X"
             "XXX"
            },
            {10,
             3, 2, 1,
             2, 2, 0,
             " X "
             "XXX"
            },
            {11,
             3, 2, 1,
             1, 2, 0,
             "X X"
             "XXX"
            },
            {12,
             3, 2, 1,
             2, 2, 0,
             " XX"
             "XXX"
            },
            {13,
             4, 1, 1,
             3, 3, 1,
             "XXXX"
            },
            {14,
             4, 2, 1,
             3, 4, 1,
             "  X "
             "XXXX"
            },
            {15,
             4, 2, 1,
             7, 4, 2,
             " XX "
             "XXXX"
            },
            {16,
             4, 2, 1,
             2, 3, 1,
             " XXX"
             "XX  "
            },
            {17,
             4, 2, 1,
             10, 3, 2,
             "   X"
             "XXXX"
            },
            {18,
             4, 2, 1,
             10, 5, 3,
             "  XX"
             "XXXX"
            },
            {19,
             4, 2, 1,
             4, 2, 0,
             " XXX"
             "XXX "
            },
            {20,
             5, 1, 1,
             7, 1, 1,
             "XXXXX"
            },
            {21,
             3, 3, 1,
             5, 4, 2,
             " X "
             "XXX"
             " X "
            },
            {22,
             3, 3, 1,
             5, 5, 2,
             "X  "
             "XXX"
             "X  "
            },
            {23,
             3, 3, 1,
             8, 6, 3,
             " XX"
             "XXX"
             "X  "
            },
            {24,
             3, 3, 1,
             3, 6, 2,
             " X "
             "XXX"
             "X X"
            },
            {25,
             3, 3, 1,
             2, 3, 0,
             "X X"
             "XXX"
             "X X"
            },
            {26,
             3, 3, 1,
             10, 4, 3,
             "  X"
             " XX"
             "XX "
            },
            {27,
             4, 3, 1,
             0, 3, 1,
             " X  "
             "XXXX"
             " X  "
            },
            {28,
             4, 3, 1,
             5, 3, 1,
             " XX "
             "XXXX"
             " XX "
            },
            {29,
             4, 3, 1,
             2, 1, 0,
             "  X "
             "XXXX"
             " X  "
            },
            {30,
             4, 3, 1,
             7, 2, 2,
             "X   "
             "XXXX"
             "X   "
            },
            {31,
             4, 3, 1,
             1, 2, 0,
             "X   "
             "XXXX"
             "   X"
            },
            {32,
             4, 2, 1,
             1, 5, 1,
             "X  X"
             "XXXX"
            },
            {33,
             5, 3, 1,
             1, 4, 1,
             "  X  "
             "XXXXX"
             "  X  "
            },

    };

    class TileSources
    {
    private:
        const std::vector<TileSource> sources_;
        const std::vector<TileSource> choosable_tile_sources_;
        const std::vector<int> choosable_tiles_;
        const std::vector<int> single_tiles_;
    public:
        static TileSources& instance()
        {
            static TileSources instance; // Guaranteed to be destroyed.
                                         // Instantiated on first use.
            return instance;
        }
    private:
        static const std::vector<TileSource> collect_sources() noexcept {
            std::vector<TileSource> result;
            for (const auto &tile : base_tiles) {
                for (int i = 0; i < tile.amount(); ++i) {
                    result.emplace_back(tile);
                }
            }
            return result;
        }

        static const std::vector<TileSource> select_choosable_sources(const std::vector<TileSource>& sources) noexcept {
            std::vector<TileSource> result;
            result.reserve(sources.size());
            for (const auto &source : sources) {
                if (source.is_start_tile()) {
                    result.emplace_back(source);
                    break;
                }
            }
            assert(result.size() == 1);
            for (const auto &source : sources) {
                if (!(source.is_bonus_tile() || source.is_start_tile())) {
                    result.emplace_back(source);
                }
            }
            return result;
        }

        static const std::vector<int> select_choosable(const std::vector<TileSource>& sources) noexcept {
            std::vector<int> result;
            result.reserve(sources.size());
            for (int i = 0; i < sources.size(); ++i) {
                if (sources[i].is_start_tile()) {
                    result.emplace_back(i);
                    break;
                }
            }
            assert(result.size() == 1);
            for (int i = 0; i < sources.size(); ++i) {
                if (!(sources[i].is_bonus_tile() || sources[i].is_start_tile())) {
                    result.emplace_back(i);
                }
            }
            return result;
        }

        static const std::vector<int> select_singles(const std::vector<TileSource>& sources) noexcept {
            std::vector<int> result;
            result.reserve(std::max(patchwork::board_singles[0], patchwork::board_singles[patchwork::board_length-1]));
            for (int i = 0; i < sources.size(); ++i) {
                if (sources[i].is_bonus_tile()) {
                    result.emplace_back(i);
                }
            }
            return result;
        }

        TileSources() :
                sources_(collect_sources()),
                choosable_tile_sources_(select_choosable_sources(sources_)),
                choosable_tiles_(select_choosable(sources_))
        {
            assert(choosable_tiles_.size() == choosable_tile_sources_.size());
        }
        friend int tile_count();
        friend int tile_count_choosable();
        friend const std::vector<int>& choosable_tiles();
        friend const std::vector<int>& single_tiles();
        friend TileSource tile(int tile);
    public:
        TileSources(TileSources const&) = delete;
        void operator=(TileSources const&)  = delete;
    };

    int tile_count() {
        return TileSources::instance().sources_.size();
    }

    int tile_count_choosable() {
        return TileSources::instance().choosable_tiles_.size();
    }

    const std::vector<int>& choosable_tiles() {
        return TileSources::instance().choosable_tiles_;
    }

    const std::vector<int>& single_tiles() {
        return TileSources::instance().single_tiles_;
    }

    TileSource tile(int tile) {
        return TileSources::instance().sources_[tile];
    }
}

#pragma clang diagnostic pop