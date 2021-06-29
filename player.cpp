#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cstdlib>
#include <ctime>
#define INF 100000
#define MINUS_INF -100000
#define DEPTH 6

using std::cout;
using std::vector;

struct Point {
    int x, y;
	Point() : Point(0, 0) {}
	Point(float x, float y) : x(x), y(y) {}
	Point(int x, int y) : x(x), y(y) {}
    bool operator==(const Point& rhs) const {
		return x == rhs.x && y == rhs.y;
	}
	bool operator!=(const Point& rhs) const {
		return !operator==(rhs);
	}
	Point operator+(const Point& rhs) const {
		return Point(x + rhs.x, y + rhs.y);
	}
	Point operator-(const Point& rhs) const {
		return Point(x - rhs.x, y - rhs.y);
	}
};
std::ostream& operator<<(std::ostream& os, Point p) {
    os << "Point(" << p.x << ", " << p.y << ")";
    return os;
}

int player;
const int SIZE = 8;
using Board = std::array<std::array<int, SIZE>, SIZE>;
Board b;
std::vector<Point> next_valid_spots;
std::ifstream fin;
std::ofstream fout;
int action_idx = 0;

class OthelloBoard {
public:
    enum SPOT_STATE {
        EMPTY = 0,
        BLACK = 1,
        WHITE = 2
    };
    static const int SIZE = 8;
    const std::array<Point, 8> directions{{
        Point(-1, -1), Point(-1, 0), Point(-1, 1),
        Point(0, -1), /*{0, 0}, */Point(0, 1),
        Point(1, -1), Point(1, 0), Point(1, 1)
    }};
    std::array<std::array<int, SIZE>, SIZE> board;
    std::vector<Point> next_valid_spots;
    std::array<int, 3> disc_count;
    int cur_player;
    bool done;
    int winner;
    int turn;
private:
    int get_next_player(int player) const {
        return 3 - player;
    }
    bool is_spot_on_board(Point p) const {
        return 0 <= p.x && p.x < SIZE && 0 <= p.y && p.y < SIZE;
    }
    int get_disc(Point p) const {
        return board[p.x][p.y];
    }
    void set_disc(Point p, int disc) {
        board[p.x][p.y] = disc;
    }
    bool is_disc_at(Point p, int disc) const {
        if (!is_spot_on_board(p))
            return false;
        if (get_disc(p) != disc)
            return false;
        return true;
    }
    bool is_spot_valid(Point center) const {
        if (get_disc(center) != EMPTY)
            return false;
        for (Point dir: directions) {
            // Move along the direction while testing.
            Point p = center + dir;
            if (!is_disc_at(p, get_next_player(cur_player)))
                continue;
            p = p + dir;
            while (is_spot_on_board(p) && get_disc(p) != EMPTY) {
                if (is_disc_at(p, cur_player))
                    return true;
                p = p + dir;
            }
        }
        return false;
    }
    void flip_discs(Point center) {
        for (Point dir: directions) {
            // Move along the direction while testing.
            Point p = center + dir;
            if (!is_disc_at(p, get_next_player(cur_player)))
                continue;
            std::vector<Point> discs({p});
            p = p + dir;
            while (is_spot_on_board(p) && get_disc(p) != EMPTY) {
                if (is_disc_at(p, cur_player)) {
                    for (Point s: discs) {
                        set_disc(s, cur_player);
                    }
                    disc_count[cur_player] += discs.size();
                    disc_count[get_next_player(cur_player)] -= discs.size();
                    break;
                }
                discs.push_back(p);
                p = p + dir;
            }
        }
    }
public:
    OthelloBoard() {
        reset();
    }
    OthelloBoard(const Board& b, int cur_player) {
        disc_count[EMPTY] = 0;
        disc_count[BLACK] = 0;
        disc_count[WHITE] = 0;
        turn = -4;
        for (int i = 0; i < SIZE; i++) {
            for (int j = 0; j < SIZE; j++) {
                board[i][j] = b[i][j];
                if(b[i][j] == EMPTY)
                    disc_count[EMPTY]++;
                else if(b[i][j] == BLACK)
                    disc_count[BLACK]++, turn++;
                else
                    disc_count[WHITE]++, turn++;
            }
        }
        this->cur_player = cur_player;
        next_valid_spots = get_valid_spots();
        done = false;
        winner = -1;
    }
    void reset() {
        for (int i = 0; i < SIZE; i++) {
            for (int j = 0; j < SIZE; j++) {
                board[i][j] = EMPTY;
            }
        }
        board[3][4] = board[4][3] = BLACK;
        board[3][3] = board[4][4] = WHITE;
        cur_player = BLACK;
        disc_count[EMPTY] = 8*8-4;
        disc_count[BLACK] = 2;
        disc_count[WHITE] = 2;
        next_valid_spots = get_valid_spots();
        done = false;
        winner = -1;
        turn = 0;
    }
    std::vector<Point> get_valid_spots() const {
        std::vector<Point> valid_spots;
        for (int i = 0; i < SIZE; i++) {
            for (int j = 0; j < SIZE; j++) {
                Point p = Point(i, j);
                if (board[i][j] != EMPTY)
                    continue;
                if (is_spot_valid(p))
                    valid_spots.push_back(p);
            }
        }
        return valid_spots;
    }
    bool put_disc(Point p) {
        if(!is_spot_valid(p)) {
            winner = get_next_player(cur_player);
            done = true;
            return false;
        }
        set_disc(p, cur_player);
        disc_count[cur_player]++;
        disc_count[EMPTY]--;
        flip_discs(p);
        // Give control to the other player.
        cur_player = get_next_player(cur_player);
        next_valid_spots = get_valid_spots();
        // Check Win
        if (next_valid_spots.size() == 0) {
            cur_player = get_next_player(cur_player);
            next_valid_spots = get_valid_spots();
            if (next_valid_spots.size() == 0) {
                // Game ends
                done = true;
                int white_discs = disc_count[WHITE];
                int black_discs = disc_count[BLACK];
                if (white_discs == black_discs) winner = EMPTY;
                else if (black_discs > white_discs) winner = BLACK;
                else winner = WHITE;
            }
        }
        return true;
    }
};

void read_board(std::ifstream& fin) {
    fin >> player;
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            fin >> b[i][j];
        }
    }
}

void read_valid_spots(std::ifstream& fin) {
    int n_valid_spots;
    fin >> n_valid_spots;
    int x, y;
    for (int i = 0; i < n_valid_spots; i++) {
        fin >> x >> y;
        next_valid_spots.push_back({x, y});
    }
}

void write_valid_spot(Point p) {
    //int n_valid_spots = next_valid_spots.size();
    //srand(time(NULL));
    
    // Choose random spot. (Not random uniform here)
    //int index = (rand() % n_valid_spots);
    //Point p = next_valid_spots[action_idx];
    // Remember to flush the output to ensure the last action is written to file.
    fout << p.x << " " << p.y << std::endl;
    fout.flush();
}

int edge(OthelloBoard &board) {
    int edge[3] = {0, 0, 0};
    
    for(int i = 1; i < 7; i++) {
        edge[board.board[i][0]]++;
    }
    for(int i = 1; i < 7; i++) {
        edge[board.board[i][7]]++;
    }
    for(int j = 1; j < 7; j++) {
        edge[board.board[0][j]]++;
    }
    for(int j = 1; j < 7; j++) {
        edge[board.board[7][j]]++;
    }
    if(edge[OthelloBoard::BLACK]+edge[OthelloBoard::WHITE] == 0)
        return 0;
    double v = ( (double) edge[OthelloBoard::BLACK]-edge[OthelloBoard::WHITE])/(edge[OthelloBoard::BLACK]+edge[OthelloBoard::WHITE]);
    return 100 * v;
}

int mobility(OthelloBoard &board) {
    int mobility[3] = {0, 0, 0};
    mobility[board.cur_player] = board.next_valid_spots.size();
    
    board.cur_player = 3 - board.cur_player;
    board.next_valid_spots = board.get_valid_spots();
    mobility[board.cur_player] = board.next_valid_spots.size();
    
    if(mobility[OthelloBoard::BLACK]+mobility[OthelloBoard::WHITE] == 0)
        return 0;
    double v = ( (double) mobility[OthelloBoard::BLACK]-mobility[OthelloBoard::WHITE])/(mobility[OthelloBoard::BLACK]+mobility[OthelloBoard::WHITE]);
    return 100 * v;
}

int potential_mobility(OthelloBoard &board) {
    int p_mobility[3] = {0, 0, 0};
    bool black_cnt, white_cnt;
    for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 8; j++) {
            if(board.board[i][j] == OthelloBoard::EMPTY) {
                black_cnt = white_cnt = true;
                if(i+1 < 8) {
                    if(white_cnt && board.board[i+1][j] == OthelloBoard::BLACK) {
                        p_mobility[OthelloBoard::WHITE]++;
                        white_cnt = false;
                    }
                    if(black_cnt && board.board[i+1][j] == OthelloBoard::WHITE) {
                        p_mobility[OthelloBoard::BLACK]++;
                        black_cnt = false;
                    }
                }
                if(i-1 >= 0) {
                    if(white_cnt && board.board[i-1][j] == OthelloBoard::BLACK) {
                        p_mobility[OthelloBoard::WHITE]++;
                        white_cnt = false;
                    }
                    if(black_cnt && board.board[i-1][j] == OthelloBoard::WHITE) {
                        p_mobility[OthelloBoard::BLACK]++;
                        black_cnt = false;
                    }
                }
                if(j+1 < 8) {
                    if(white_cnt && board.board[i][j+1] == OthelloBoard::BLACK) {
                        p_mobility[OthelloBoard::WHITE]++;
                        white_cnt = false;
                    }
                    if(black_cnt && board.board[i][j+1] == OthelloBoard::WHITE) {
                        p_mobility[OthelloBoard::BLACK]++;
                        black_cnt = false;
                    }
                }
                if(j-1 >= 0) {
                    if(white_cnt && board.board[i][j-1] == OthelloBoard::BLACK) {
                        p_mobility[OthelloBoard::WHITE]++;
                        white_cnt = false;
                    }
                    if(black_cnt && board.board[i][j-1] == OthelloBoard::WHITE) {
                        p_mobility[OthelloBoard::BLACK]++;
                        black_cnt = false;
                    }
                }
            }
        }
    }
    if(p_mobility[OthelloBoard::BLACK]+p_mobility[OthelloBoard::WHITE] == 0)
        return 0;
    double v = ( (double) p_mobility[OthelloBoard::BLACK]-p_mobility[OthelloBoard::WHITE]) / (p_mobility[OthelloBoard::BLACK]+p_mobility[OthelloBoard::WHITE]);
    return 100 * v;
    
}

int stability(OthelloBoard &board) {
    int stability[3] = {0, 0, 0};
    int upper_left = board.board[0][0];
    int upper_right = board.board[0][7];
    int lower_left = board.board[7][0];
    int lower_right = board.board[7][7];
    bool calc_reverse[4] = {};
    
    // they are counted once more
    stability[lower_left]--, stability[lower_right]--;
    stability[upper_left]--, stability[upper_right]--;
 
    for(int j = 0; j < 8; j++) { // upper left to upper right
        if(board.board[0][j] == upper_left) {
            stability[upper_left]++;
        }
        else {
            calc_reverse[0] = true;
            break;
        }
    }
    for(int j = 7; calc_reverse[0] && j >= 0; j--) { // upper right to upper left
        if(board.board[0][j] == upper_right) {
            stability[upper_right]++;
        }
        else {
            break;
        }
    }
    for(int i = 0; i < 8; i++) { // upper left to lower left
        if(board.board[i][0] == upper_left) {
            stability[upper_left]++;
        }
        else {
            calc_reverse[1] = true;
            break;
        }
    }
    for(int i = 7; calc_reverse[1] && i >= 0; i--) { // lower left to upper left
        if(board.board[i][0] == lower_left) {
            stability[lower_left]++;
        }
        else {
            break;
        }
    }
    for(int j = 0; j < 8; j++) { // lower left to lower right
        if(board.board[7][j] == lower_left) {
            stability[lower_left]++;
        }
        else {
            calc_reverse[2] = true;
            break;
        }
    }
    for(int j = 7; calc_reverse[2] && j >= 0; j--) { // lower right to lower left;
        if(board.board[7][j] == lower_right) {
            stability[lower_right]++;
        }
        else {
            break;
        }
    }
    for(int i = 0; i < 8; i++) { // upper right to lower right
        if(board.board[i][7] == upper_right) {
            stability[upper_right]++;
        }
        else {
            calc_reverse[3] = true;
            break;
        }
    }
    for(int i = 7; calc_reverse[3] && i >= 0; i--) { // lower right to upper right
        if(board.board[i][7] == lower_right) {
            stability[lower_right]++;
        }
        else {
            break;
        }
    }
    if(stability[OthelloBoard::BLACK]+stability[OthelloBoard::WHITE] == 0)
        return 0;
    double v = ( (double) stability[OthelloBoard::BLACK]-stability[OthelloBoard::WHITE]) / (stability[OthelloBoard::BLACK]+stability[OthelloBoard::WHITE]);
    return 100 * v;
}

int corner(OthelloBoard &board) {
    int corner[3] = {0, 0, 0};
    corner[board.board[0][0]]++, corner[board.board[0][7]]++;
    corner[board.board[7][0]]++, corner[board.board[7][7]]++;
    if(corner[OthelloBoard::BLACK]+corner[OthelloBoard::WHITE] == 0)
        return 0;
    double v = ( (double) corner[OthelloBoard::BLACK]-corner[OthelloBoard::WHITE]) / (corner[OthelloBoard::BLACK]+corner[OthelloBoard::WHITE]);
    return 100 * v;
}

int coin_parity(OthelloBoard &board) {
    int coin_parity[3] = {0, board.disc_count[1], board.disc_count[2]};
    if(coin_parity[OthelloBoard::BLACK]+coin_parity[OthelloBoard::WHITE] == 0)
        return 0;
    double v = ( (double) coin_parity[OthelloBoard::BLACK]-coin_parity[OthelloBoard::WHITE]) / (coin_parity[OthelloBoard::BLACK]+coin_parity[OthelloBoard::WHITE]);
    return 100 * v;
}

int evaluate(OthelloBoard &board) {
    int v;
    if(board.turn < 15) {
        v = 25 * corner(board) +
            100 * mobility(board) +
            100 * potential_mobility(board) +
            25 * coin_parity(board) +
            10 * edge(board) +
            25 * stability(board);
    }
    else if(board.turn < 55) {
        v = 100 * corner(board) + 
            25 * mobility(board) +
            25 * potential_mobility(board) +
            20 * edge(board) +
            //25 * coin_parity(board) + 
            100 * stability(board);
    }
    else {
        v = 100 * corner(board) + 
            //5 * mobility(board) +
            200 * coin_parity(board) +
            100 * stability(board);
    }

    if(player == OthelloBoard::BLACK)
        return v;
    else
        return -v;
}

int minimax(OthelloBoard &board, int depth, int alpha, int beta, int cur_player) {
    if(!depth || board.next_valid_spots.empty()) {
        //cout << "\t\tscore: " << evaluate(board) << "\n";
        return evaluate(board);
    }
    int v, tmp_score;
    //cout << "calling minimax with depth " << depth << "\n";
    if(cur_player == player) { // max node
        //cout << "\tthis is max node\n\ttraversing children:\n";
        v = MINUS_INF;
        for(size_t i = 0; i < board.next_valid_spots.size(); i++) {
            auto p = board.next_valid_spots[i];
            //cout << "\t\tnext step: " << p << "\n";
            auto next_board = board;
            next_board.put_disc(p);
            tmp_score = minimax(next_board, depth-1, alpha, beta, 3-cur_player);
            if(tmp_score > v) {
                v = tmp_score;
                if(depth == DEPTH) {
                    //cout << "back to original depth and update: " << board.next_valid_spots[i] << "\n";
                    //cout << "score: " << tmp_score << "\n";
                    action_idx = i;
                    write_valid_spot(board.next_valid_spots[i]);
                }
            }
            else if(depth == DEPTH) {
                //cout << "back to original depth\n";
                //cout << board.next_valid_spots[i] << "is not good enough\n";
                //cout << "score: " << tmp_score << "\n";

            }
            alpha = std::max(v, alpha);
            if(beta <= alpha) {
                //cout << "alpha pruning\n";
                break;
            }
        }
    }
    else { // min node
        v = INF;
        //cout << "\tthis is min node\n\ttraversing children:\n";
        for(size_t i = 0; i < board.next_valid_spots.size(); i++) {
            auto p = board.next_valid_spots[i];
            //cout << "\t\tnext step: " << p << "\n";
            auto next_board = board;
            next_board.put_disc(p);
            v = std::min(v, minimax(next_board, depth-1, alpha, beta, 3-cur_player));
            beta = std::min(v, beta);
            if(beta <= alpha) {
                //cout << "beta pruning\n";
                break;
            }
        }
    }
    return v;
}

int main(int, char** argv) {
    fin = std::ifstream(argv[1]);
    fout = std::ofstream(argv[2]);
    read_board(fin);
    read_valid_spots(fin);

    OthelloBoard board(b, player);
    minimax(board, DEPTH, MINUS_INF, INF, player);
    //cout << evaluate(board);
    fin.close();
    fout.close();
    return 0;
}
