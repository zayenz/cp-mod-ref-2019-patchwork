//
// Created by Mikael Zayenz Lagerkvist
//

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCUnusedGlobalDeclarationInspection"
#ifndef PATCHWORK_TILES_H
#define PATCHWORK_TILES_H

#include <vector>
#include <gecode/minimodel.hh>

namespace patchwork {

    /** \brief Specification of one tile
     *
     * This structure can be used to specify a tile with specified width
     * and height, number of such tiles (all with unique values), and a
     * char-array tile showing the tile in row-major order.x
     *
     * \relates PatchworkModel
     */
    class Tile {
        const int width;  ///< Width of tile
        const int height; ///< Height of tile
        const std::vector<bool> marks; ///< The row-major marks of the tile in a minimum bounding box
    public:
        Tile(int width, int height, std::vector<bool> marks);

        const Gecode::REG make_placement_expression(const Gecode::REG& mark, const Gecode::REG& empty_mark, int board_width, int board_height) const;
        const Gecode::REG make_usage_expression() const;

        inline bool operator()(int x, int y) const {
            return at(x, y);
        }

        /**
         *
         * @param x The x-coordinate (the column)
         * @param y The y-coordinate (the row)
         * @return The value at (x, y)
         */
        inline bool at(int x, int y) const {
            assert(0 <= x && x < width);
            assert(0 <= y && y < height);
            return marks[y * width + x];
        }

        inline bool operator==(const Tile &rhs) const noexcept {
            if (width != rhs.width || height != rhs.height) {
                return false;
            }
            for (int i = 0; i < marks.size(); ++i) {
                if (marks[i] != rhs.marks[i]) {
                    return false;
                }
            }
            return true;
        }

        inline bool operator!=(const Tile &rhs) const noexcept {
            return !(*this == rhs);
        }

    private:
        const std::vector<int> compute_column_usage() const noexcept;
        const std::vector<int> compute_row_usage() const noexcept;
    };

    /** \brief Specification of one tile
     *
     * This structure can be used to specify a tile with specified width
     * and height, number of such tiles (all with unique values), and a
     * char-array tile showing the tile in row-major order.
     *
     * \relates PatchworkModel
     */
    class TileSource {
        const int id_; ///< The id for this source
        const int amount_; ///< Number of copies of the tile
        const int area_; ///< The area occupied by the tile
        const int button_cost_; ///< The button cost for the tile
        const int time_cost_; ///< The time cost for the tile
        const int buttons_; ///< The number of buttons income for the tile
        const std::vector<Tile> alternatives_; ///< The tile alternatives (id, rot90, flip vertical, ...)
        const Gecode::REG placement_expression_; ///< A placement expression for placing on a 9x9 grid of Boolean variables with 1 eol column
        const Gecode::REG usage_expression_; ///< A usage expression for 9 columns followed by an empty end-of-line column and then 9 rows
    public:
        /**
         *
         * @param width The width of the tile pattern
         * @param height The height of the tile pattern
         * @param amount The number of copies of the tile
         * @param button_cost The button cost for the tile
         * @param time_cost The time cost for the tile
         * @param buttons The number of buttons of income for the tile
         * @param tile_pattern A row-major tile pattern
         */
        TileSource(int id, int width, int height, int amount, int button_cost,
                   int time_cost, int buttons, const char *tile_pattern);

        /**
         *
         * @return A placement expression for the tile on a 9x9 grid of Boolean variables with 1 eol column.
         */
        const Gecode::REG as_placement_expression() const;

        /**
         *
         * @return A usage expression for the tile for the columns followed by the rows
         */
        const Gecode::REG usage_expression() const;

        const int id() const;

        const int amount() const;

        const int area() const;

        const int button_cost() const;

        const int time_cost() const;

        const int buttons() const;

        /**
         *
         * @return True iff this tile is the tile to start after
         */
        const bool is_start_tile() const;

        /**
         *
         * @return True iff the tile is a bonus-tile from the board.
         */
        const bool is_bonus_tile() const;

        /**
         * Give an evaluation fo the value of this tile when considering buying it at a certain square.
         *
         * @param square The square in which to evaluate
         * @return The value of this tile type
         */
        const double value(int square) const;

    private:
        /**
         *
         * @param width The width of the tile pattern
         * @param height The height of the tile pattern
         * @param tile_pattern A row-major tile pattern
         * @return A vector of all the uniquely transformed tiles
         */
        static std::vector<Tile> make_unique_tiles(int width, int height, const char *tile_pattern);

        /**
         *
         * @param width The width of the tile pattern
         * @param height The height of the tile pattern
         * @param tile_pattern A row-major tile pattern
         * @return The area occupied by the tile
         */
        static const int count_area(int width, int height, const char *pattern);

        /**
         * Create the expression for placing this tile, with the alternative index as the first 8 variables.
         *
         * @param alternatives The alternatives to make the expression over.
         * @return A placement expression on a 9x9 grid with 1 eol column, preceeded by the alternative index
         */
        static const Gecode::REG make_placement_expression(const std::vector<Tile>& alternatives);

        /**
         *
         * @param alternatives The alnernatives to make the expression over
         * @return An usage expression on 9 column usages, 1 separator, and 9 row usages, preceeded by the alternative index
         */
        static const Gecode::REG make_usage_expression(const std::vector<Tile>& alternatives);

        /**
         * The alternative index expression is an expression over Boolean variables indicating which, alternative among the 8 possible it is.
         *
         * An expression consists of a zero, alternative-1 zeroes, a 1, and then 8-alternative zeroes. The leading zero is due to the choice
         * to make the alternative indexes 1-based in the expression, to make non-placement (that is, no alternative used) a possibility.
         *
         * @param i The alternative index (between 0 and 7)
         * @return The alternative index expressions
         */
        static const Gecode::REG make_alternative_index(int alternative);
    };

    /**
     * All the tiles in Patchwork
     */
    extern const std::vector<TileSource> base_tiles;

    /**
     *
     * @return The number of unique tiles
     */
    int tile_count();

    /**
     *
     * @return The number of unique tiles that are choosable (that is, not including the 1x1 tiles)
     */
    int tile_count_choosable();

    /**
     * Guarantees that that the tile in postion 0 is the start tile, and the rest are all tiles to choose from
     * and not the 1x1 bonus tiles.
     *
     * @return The unique tiles that are choosable (that is, not including the 1x1 tiles)
     */
    const std::vector<int>& choosable_tiles();

    /**
     * Only the bonus tiles.
     *
     * @return The unique bonus/single tiles.
     */
    const std::vector<int>& single_tiles();

    /**
     *
     * @param tile The tile index to get. Must be between 0 and tile_count()
     * @return The source for tile number \a tile
     */
    TileSource tile(int tile);
}

#endif //PATCHWORK_TILES_H

#pragma clang diagnostic pop