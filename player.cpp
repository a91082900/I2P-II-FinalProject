#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <unordered_map>
#include <cstdlib>
#include <ctime>
#define INF 100000
#define MINUS_INF -100000
#define INIT_DEPTH 3
#define MAX_DEPTH 6

using std::cout;
using std::vector;

unsigned long long hash_val[8][8][3] = {
    0x7c9befc91963be48ULL, 9852353736228043306ULL, 6958017284889570045ULL, 
    0xbb7830a71d4e0246ULL, 11589589853262412901ULL, 10315806072214505521ULL, 
    0x830026e163c9850aULL, 10818023254006820513ULL, 17572290972824691223ULL, 
    0x9a95e4781ff7db57ULL, 458407484403071964ULL, 16652370706459499591ULL, 
    0x9d5a98b7fbb65b45ULL, 16500692775934581219ULL, 5790415578591874884ULL, 
    0x506aeb105c71cd2bULL, 4759439104965278017ULL, 13128712239172345819ULL, 
    0x777766b514ed3c41ULL, 2295549573310786278ULL, 7603549769048207158ULL, 
    0x48a030bf63970903ULL, 10225568473374670931ULL, 9059832052614340196ULL, 
    0x4042c85f4428b0fULL, 11328474593544583058ULL, 11172701540615104828ULL, 
    0x20d70d69bfaa0f3fULL, 6185611776103876737ULL, 5734999815514311146ULL, 
    0xb82c1fb03e63961fULL, 12488028590855641894ULL, 11635725213107818965ULL, 
    0xc7c6a93597b11417ULL, 10492333263082880106ULL, 5507212678140400363ULL, 
    0xc0a1729714470adcULL, 15112231158906992524ULL, 8804976461615772751ULL, 
    0x890b4eb18b08519ULL, 15281768490874243036ULL, 16133627091960451837ULL, 
    0xf44b155d574797cdULL, 9142155325104945848ULL, 8495916523671078207ULL, 
    0xbbd31ef2d9a307bdULL, 15302987277423351893ULL, 15305266803211037826ULL, 
    0xb2e96e8f516a0babULL, 15309201425878274502ULL, 16708669793962872707ULL, 
    0xd7488a69da3cd8b7ULL, 6076201160760333852ULL, 14901759792555822931ULL, 
    0xc1f0dba8059a7a56ULL, 14121309066462108791ULL, 5329416841205627558ULL, 
    0xbc32708311fadd88ULL, 5343944193514868508ULL, 16121597139585075523ULL, 
    0xca09872b668741bcULL, 5866279830796149579ULL, 3554437239745204579ULL, 
    0x5d2c96878abf8882ULL, 7899073867915714743ULL, 13789264105200029724ULL, 
    0xe56e858bbeae70aaULL, 3996389219710758659ULL, 8215150355406308188ULL, 
    0x20952230581da701ULL, 14646177331817167761ULL, 12158548635227981911ULL, 
    0x2e490033324827ecULL, 826179905950230131ULL, 16089344200505845216ULL, 
    0xdcd7ecf356e157c8ULL, 16878253304889185728ULL, 9637422128555754172ULL, 
    0x1338458387e7c7d0ULL, 10563120243578066817ULL, 5869084023780218714ULL, 
    0x11ea8c8e5a37ed1fULL, 1991643539980122473ULL, 10557767269249362051ULL, 
    0x3843d40dad97ba82ULL, 17143129812335491161ULL, 6431067239510716949ULL, 
    0x78d227f1dd32207fULL, 7851314334235298928ULL, 5477029574442150210ULL, 
    0xa7e4b6d1fb6069b9ULL, 14553877516163379115ULL, 2363364968091343983ULL, 
    0xbf99bfb93ee6d322ULL, 5677840183623434615ULL, 66886679689720648ULL, 
    0xbc6ac69c1c2a91faULL, 5444231061531287163ULL, 9695006787525442380ULL, 
    0xe6f7f8206a10f303ULL, 18305080048128855577ULL, 5645267922000355504ULL, 
    0x5c43631d2acf704aULL, 12929928375785092746ULL, 4377898355178459819ULL, 
    0x862518568be1cbe7ULL, 16187017735932326025ULL, 13030818728300357728ULL, 
    0x25f2549d4845d4ecULL, 10959725368111771313ULL, 3691567413155767177ULL, 
    0xaba2e567265b7fa0ULL, 1800031120175556009ULL, 18405348092835211913ULL, 
    0xe905c5a04195b2f4ULL, 12730539848014610468ULL, 8333056188286962192ULL, 
    0x3a5d95d883fc5e6bULL, 16176519059140742957ULL, 2455479682650861727ULL, 
    0xfb05efb016f9015cULL, 12633670182190845216ULL, 14725343082164666053ULL, 
    0xfef84769bb931783ULL, 256853796529907756ULL, 13588571900527738334ULL, 
    0xe1a99616e44f00dfULL, 135232709750439629ULL, 8474155065250475884ULL, 
    0xd379e49947e37347ULL, 5372159243758834570ULL, 17015132724908581253ULL, 
    0xc07b18435ffea521ULL, 6531202999527767822ULL, 8336198634133900512ULL, 
    0x99c268d112c6b07ULL, 7886496841473319076ULL, 17292526927850174375ULL, 
    0x261603e80593c1a8ULL, 9739106225599917966ULL, 2780501810267444869ULL, 
    0xf6cf2a93feab4ad7ULL, 10685845474162025124ULL, 15340955385287392480ULL, 
    0x3e6bdd57d487fd1cULL, 13978964174790994831ULL, 8644273313796362959ULL, 
    0x1e140ed508879081ULL, 16022747951384283742ULL, 15723834684104944235ULL, 
    0xe18e8d21baf48df7ULL, 6197291079413624081ULL, 17435474805534094702ULL, 
    0x117d57d8d491eccfULL, 6227608893576468309ULL, 13770141428937724979ULL, 
    0xeaba688a54d74200ULL, 3184610583047198441ULL, 9317370869037172179ULL, 
    0x5f834b24d6ec4b88ULL, 9328283286673496172ULL, 4956093264878840044ULL, 
    0x49e30da0a2f234faULL, 16111719831761529677ULL, 5305188443398841454ULL, 
    0x38f3364fda0bfeaeULL, 14988040926854518015ULL, 2138526800182427739ULL, 
    0x7693c3467e6d8dfcULL, 5528415974937514314ULL, 16951075450268602148ULL, 
    0xe80d0d902bce809bULL, 8845130567742619749ULL, 17051509938282636425ULL, 
    0x6299399d094c4a0fULL, 11166558137496532999ULL, 6497153536099740346ULL, 
    0x1fb175b252fad0c8ULL, 16951643777001253420ULL, 13346740626772796121ULL, 
    0xb2e87247eeaa48fbULL, 12977341909087805581ULL, 183491818910736501ULL, 
    0x4e993f0ee91b0b4bULL, 4215405735574049812ULL, 446488715112696539ULL, 
    0xcb103d2079306433ULL, 10506199991675698391ULL, 5538179496161336337ULL, 
    0x6224d382f3b9affeULL, 6478065130795735537ULL, 10199748655594304210ULL
};

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
int now_depth;
Point write_point;

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
    unsigned long long hash = 0;
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
        hash ^= hash_val[p.x][p.y][board[p.x][p.y]];
        hash ^= hash_val[p.x][p.y][disc];
        board[p.x][p.y] = disc;
    }
    void set_disc(int x, int y, int disc) {
        Point p(x, y);
        set_disc(p, disc);
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
        hash = 0;
        for (int i = 0; i < SIZE; i++) {
            for (int j = 0; j < SIZE; j++) {
                board[i][j] = b[i][j];
                hash ^= hash_val[i][j][b[i][j]];
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
        hash = 0;
        for (int i = 0; i < SIZE; i++) {
            for (int j = 0; j < SIZE; j++) {
                board[i][j] = EMPTY;
                hash ^= hash_val[i][j][EMPTY];
            }
        }
        set_disc(3, 4, BLACK);
        set_disc(4, 3, BLACK);
        set_disc(3, 3, WHITE);
        set_disc(4, 4, BLACK);

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

void write_valid_spot() {
    //int n_valid_spots = next_valid_spots.size();
    //srand(time(NULL));
    // Choose random spot. (Not random uniform here)
    //int index = (rand() % n_valid_spots);
    //Point p = next_valid_spots[action_idx];
    // Remember to flush the output to ensure the last action is written to file.
    fout << write_point.x << " " << write_point.y << std::endl;
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
            75 * mobility(board) +
            75 * potential_mobility(board) +
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
                if(depth == now_depth) {
                    //cout << "back to original depth and update: " << board.next_valid_spots[i] << "\n";
                    //cout << "score: " << tmp_score << "\n";
                    write_point = board.next_valid_spots[i];
                }
            }
            else if(depth == now_depth) {
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

struct TTNode {
    int depth;
    int value;
    enum FLAG_TYPE {
        UPPERBOUND, LOWERBOUND, EXACT
    };
    int flag;
};

std::unordered_map<unsigned long long, TTNode> TT;
int minimax_with_tt(OthelloBoard &board, int depth, int alpha, int beta, int cur_player) {
    int orig_a = alpha, orig_b = beta;
    auto entry = TT.find(board.hash);
    TTNode tt_entry;
    if(entry != TT.end()) {
        tt_entry = entry->second;
        if(tt_entry.depth >= depth) {
            if(tt_entry.flag == TTNode::EXACT)
                return tt_entry.value;
            else if(tt_entry.flag == TTNode::UPPERBOUND) {
                beta = std::min(beta, tt_entry.value);
            }
            else {
                alpha = std::max(alpha, tt_entry.value);
            }
        }
        if (beta <= alpha)
            return tt_entry.value;
    }
    //cout << board.hash << "\n";
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
            tmp_score = minimax_with_tt(next_board, depth-1, alpha, beta, 3-cur_player);
            if(tmp_score > v) {
                v = tmp_score;
                if(depth == now_depth) {
                    //cout << "back to original depth and update: " << board.next_valid_spots[i] << "\n";
                    //cout << "score: " << tmp_score << "\n";
                    write_point = board.next_valid_spots[i];
                }
            }
            else if(depth == now_depth) {
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
            v = std::min(v, minimax_with_tt(next_board, depth-1, alpha, beta, 3-cur_player));
            beta = std::min(v, beta);
            if(beta <= alpha) {
                //cout << "beta pruning\n";
                break;
            }
        }
    }

    tt_entry.value = v;
    if(v <= orig_a) {
        tt_entry.flag = TTNode::UPPERBOUND;
    }
    else if(v >= orig_b) {
        tt_entry.flag = TTNode::LOWERBOUND;
    }
    else
        tt_entry.flag = TTNode::EXACT;
    tt_entry.depth = depth;	
    TT[board.hash] = tt_entry;
    return v;
}

int main(int, char** argv) {
    fin = std::ifstream(argv[1]);
    fout = std::ofstream(argv[2]);
    read_board(fin);
    read_valid_spots(fin);

    OthelloBoard board(b, player);
    now_depth = INIT_DEPTH;
    while(now_depth <= MAX_DEPTH) {
        minimax_with_tt(board, now_depth, MINUS_INF, INF, player);
        write_valid_spot();
        //cout << "finish depth " << now_depth << "\n";
        now_depth++;
    }
    fin.close();
    fout.close();
    return 0;
}
